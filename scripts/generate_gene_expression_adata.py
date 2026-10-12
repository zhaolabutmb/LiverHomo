import argparse
import os
import sys
import tempfile
from pathlib import Path

import anndata as ad
import h5py

PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.chdir(PROJECT_ROOT)
sys.path.insert(0, str(PROJECT_ROOT))

from env_config import load_env  # noqa: E402

load_env()

OBS_COLUMNS = ("cell_type2", "diagnosis", "total_counts")


def _set_column_order(group: h5py.Group, columns: tuple[str, ...]) -> None:
    group.attrs["column-order"] = list(columns)


def build_gene_expression_adata(input_path: Path, output_path: Path, force: bool = False) -> None:
    if input_path.resolve() == output_path.resolve():
        raise ValueError("Input and output AnnData paths must be different.")
    if not input_path.is_file():
        raise FileNotFoundError(f"Source AnnData not found: {input_path}")
    if output_path.exists() and not force:
        raise FileExistsError(f"Output already exists: {output_path} (pass --force to replace it)")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    temp_path = None
    try:
        with tempfile.NamedTemporaryFile(
            prefix=f".{output_path.name}.",
            suffix=".tmp",
            dir=output_path.parent,
            delete=False,
        ) as temp_file:
            temp_path = Path(temp_file.name)

        with h5py.File(input_path, "r") as source, h5py.File(temp_path, "w") as target:
            if "raw/X" not in source or "raw/var" not in source:
                raise ValueError("Source AnnData must include the full-gene raw expression matrix.")
            if "layers/lognorm" not in source:
                raise ValueError("Source AnnData must include the normalized log-expression layer.")
            if "gene_name" not in source["var"] or "gene_name" not in source["raw/var"]:
                raise ValueError("Source AnnData must include gene_name metadata for main and raw genes.")
            missing_obs = [column for column in OBS_COLUMNS if f"obs/{column}" not in source]
            if missing_obs:
                raise ValueError(f"Source AnnData is missing required observation fields: {', '.join(missing_obs)}")
            for key, value in source.attrs.items():
                target.attrs[key] = value

            source.copy("X", target, name="X")
            source.copy("var", target, name="var")
            for name in list(target["var"].keys()):
                if name not in {"_index", "gene_name"}:
                    del target["var"][name]
            _set_column_order(target["var"], ("gene_name",))

            layers = target.create_group("layers")
            layers.attrs.update(source["layers"].attrs)
            source.copy("layers/lognorm", layers, name="lognorm")

            source.copy("raw", target, name="raw")
            for name in list(target["raw/var"].keys()):
                if name not in {"_index", "gene_name"}:
                    del target[f"raw/var/{name}"]
            _set_column_order(target["raw/var"], ("gene_name",))
            del target["raw/varm"]

            source.copy("obs", target, name="obs")
            for name in list(target["obs"].keys()):
                if name != "_index" and name not in OBS_COLUMNS:
                    del target["obs"][name]
            _set_column_order(target["obs"], OBS_COLUMNS)

        with h5py.File(temp_path, "r") as check:
            expected_groups = {"X", "obs", "var", "layers", "raw"}
            if set(check.keys()) != expected_groups:
                raise ValueError(f"Generated AnnData contains unexpected groups: {sorted(check.keys())}")
            if list(check["X"].attrs["shape"]) != [
                check["obs/_index"].shape[0],
                check["var/_index"].shape[0],
            ]:
                raise ValueError("Expression matrix dimensions do not match observation and gene metadata.")

        compact = ad.read_h5ad(temp_path, backed="r")
        try:
            if compact.raw is None or set(compact.layers) != {"lognorm"} or compact.obsm:
                raise ValueError("Generated AnnData contains data outside the compact expression view.")
            if not set(OBS_COLUMNS).issubset(compact.obs.columns) or "gene_name" not in compact.var.columns:
                raise ValueError("Generated AnnData is missing required gene-expression metadata.")
        finally:
            compact.file.close()

        if output_path.exists() and not force:
            raise FileExistsError(f"Output already exists: {output_path} (pass --force to replace it)")
        os.chmod(temp_path, input_path.stat().st_mode & 0o666)
        os.replace(temp_path, output_path)
        temp_path = None
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Create an AnnData file containing full-gene expression and required gene-viewer metadata."
    )
    parser.add_argument("--input", type=Path, default=os.environ.get("DATA_PATH"))
    parser.add_argument("--output", type=Path, default=os.environ.get("GENE_EXPRESSION_DATA_PATH"))
    parser.add_argument("--force", action="store_true", help="Replace an existing output file.")
    args = parser.parse_args()
    if args.input is None:
        parser.error("Pass --input or configure DATA_PATH.")
    if args.output is None:
        args.output = args.input.parent / "adata_gene_expression.h5ad"
    build_gene_expression_adata(args.input, args.output, force=args.force)
    print(f"Generated compact gene-expression AnnData: {args.output}")


if __name__ == "__main__":
    main()
