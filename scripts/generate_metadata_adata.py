import argparse
import os
import sys
import tempfile
from pathlib import Path

import anndata as ad
import h5py
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]
os.chdir(PROJECT_ROOT)
sys.path.insert(0, str(PROJECT_ROOT))

from env_config import load_env  # noqa: E402

load_env()


def build_metadata_adata(input_path: Path, output_path: Path, force: bool = False) -> None:
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

        empty = ad.AnnData(
            obs=pd.DataFrame(index=pd.Index([], dtype=str)),
            var=pd.DataFrame(index=pd.Index([], dtype=str)),
        )
        empty.write_h5ad(temp_path)

        with h5py.File(input_path, "r") as source:
            if "obs/_index" not in source:
                raise ValueError("Source AnnData must include observation metadata.")
            expected_rows = source["obs/_index"].shape[0]
            with h5py.File(temp_path, "r+") as target:
                del target["obs"]
                source.copy("obs", target, name="obs")

        metadata = ad.read_h5ad(temp_path, backed="r")
        try:
            if metadata.n_obs != expected_rows or metadata.n_vars != 0 or "X" in metadata.file:
                raise ValueError("Generated metadata AnnData must contain observations only.")
            if metadata.obs.shape[1] == 0 and expected_rows:
                raise ValueError("Generated metadata AnnData is missing its observation columns.")
        finally:
            metadata.file.close()

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
        description="Create an observation-only AnnData file with all cell metadata rows."
    )
    parser.add_argument("--input", type=Path, default=os.environ.get("DATA_PATH"))
    parser.add_argument("--output", type=Path, default=os.environ.get("METADATA_PATH"))
    parser.add_argument("--force", action="store_true", help="Replace an existing output file.")
    args = parser.parse_args()
    if args.input is None:
        parser.error("Pass --input or configure DATA_PATH.")
    if args.output is None:
        args.output = args.input.parent / "adata_metadata.h5ad"
    build_metadata_adata(args.input, args.output, force=args.force)
    print(f"Generated compact metadata AnnData: {args.output}")


if __name__ == "__main__":
    main()
