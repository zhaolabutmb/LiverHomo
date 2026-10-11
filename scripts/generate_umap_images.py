import argparse
import json
import os
import sys
import tempfile
from pathlib import Path

import h5py
import matplotlib
import numpy as np
import pandas as pd
from anndata.io import read_elem
from dotenv import load_dotenv

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from umap_config import UMAP_METADATA_COLUMNS, UMAP_STATIC_IMAGE_FILES


def write_json_atomic(path: Path, payload: dict) -> None:
    temp_path = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=path.parent, delete=False) as handle:
            temp_path = Path(handle.name)
            json.dump(payload, handle, indent=2)
            handle.write("\n")
        os.chmod(temp_path, 0o644)
        os.replace(temp_path, path)
    finally:
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)


def render_umap(coords: np.ndarray, values: pd.Series, column: str, path: Path) -> dict:
    fig, ax = plt.subplots(figsize=(12, 8))
    temp_path = None
    try:
        numeric = pd.api.types.is_numeric_dtype(values.dtype) and not pd.api.types.is_bool_dtype(values.dtype)
        if numeric:
            colors = values.to_numpy(dtype=float, na_value=np.nan)
            colors[~np.isfinite(colors)] = np.nan
            cmap = plt.get_cmap("viridis").with_extremes(bad="#bdbdbd")
            points = ax.scatter(
                coords[:, 0], coords[:, 1], c=np.ma.masked_invalid(colors),
                cmap=cmap, s=0.5, linewidths=0, plotnonfinite=True,
            )
            fig.colorbar(points, ax=ax, label=column, fraction=0.035, pad=0.03)
            categories = None
        else:
            labels = values.astype("string").fillna("Unknown")
            categories = sorted(labels.unique().tolist())
            codes = pd.Categorical(labels, categories=categories).codes
            palette = plt.get_cmap("tab20" if len(categories) <= 20 else "gist_rainbow")
            colors = palette(np.linspace(0, 1, len(categories)))
            ax.scatter(coords[:, 0], coords[:, 1], c=colors[codes], s=0.5, linewidths=0)
            handles = [
                Line2D([], [], marker="o", linestyle="", color=color, label=category, markersize=6)
                for category, color in zip(categories, colors)
            ]
            ax.legend(handles=handles, title=column, loc="upper left", bbox_to_anchor=(1.02, 1), fontsize=8)
        ax.set_title(f"UMAP by {column}\nAll cell types and diseases ({len(coords):,} cells)")
        ax.set_xlabel("UMAP 1")
        ax.set_ylabel("UMAP 2")
        ax.set_aspect("equal", adjustable="box")
        ax.grid(False)
        fig.tight_layout()
        path.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(dir=path.parent, suffix=".png", delete=False) as handle:
            temp_path = Path(handle.name)
        fig.savefig(temp_path, format="png", dpi=160)
        os.chmod(temp_path, 0o644)
        os.replace(temp_path, path)
        return {
            "column": column,
            "cells_plotted": len(coords),
            "color_type": "numeric" if numeric else "categorical",
            "missing_values": int(values.isna().sum()),
            "categories": categories,
        }
    finally:
        plt.close(fig)
        if temp_path is not None:
            temp_path.unlink(missing_ok=True)


def generate_umap_images(input_path: Path, figure_root: Path, catalog_path: Path) -> dict:
    with catalog_path.open(encoding="utf-8") as handle:
        catalog = json.load(handle)
    plots = []
    with h5py.File(input_path, "r") as source:
        if "obsm/X_umap" not in source:
            raise ValueError("Source AnnData does not contain saved UMAP coordinates (obsm/X_umap).")
        missing = [column for column in UMAP_METADATA_COLUMNS if column not in source["obs"]]
        if missing:
            raise ValueError(f"Source AnnData is missing UMAP metadata columns: {', '.join(missing)}")
        coords = np.asarray(source["obsm/X_umap"])
        if coords.ndim != 2 or coords.shape[1] != 2 or len(coords) == 0 or not np.isfinite(coords).all():
            raise ValueError("UMAP coordinates must be a non-empty, finite, two-column matrix.")
        for column in UMAP_METADATA_COLUMNS:
            values = pd.Series(read_elem(source["obs"][column]))
            if len(values) != len(coords):
                raise ValueError(f"Metadata length does not match UMAP coordinates for {column}.")
            path = figure_root / UMAP_STATIC_IMAGE_FILES[column]
            record = render_umap(coords, values, column, path)
            record["file"] = UMAP_STATIC_IMAGE_FILES[column]
            plots.append(record)
            print(f"Generated {column}: {len(coords):,} cells -> {path}", flush=True)
    manifest = {"cells": len(coords), "subsampled": False, "filters": {}, "plots": plots}
    write_json_atomic(figure_root / "umap" / "manifest.json", manifest)
    catalog["umap_columns"] = list(UMAP_METADATA_COLUMNS)
    write_json_atomic(catalog_path, catalog)
    return manifest


def main() -> None:
    load_dotenv(PROJECT_ROOT / ".env")
    parser = argparse.ArgumentParser(description="Generate static, full-dataset UMAPs for every UMAP_METADATA_COLUMNS entry.")
    parser.add_argument("--input", type=Path, default=os.environ.get("DATA_PATH"))
    parser.add_argument("--figures", type=Path, default=PROJECT_ROOT / "figs")
    parser.add_argument("--catalog", type=Path, default=PROJECT_ROOT / "static" / "filter_catalog.json")
    args = parser.parse_args()
    if args.input is None:
        parser.error("Pass --input or configure DATA_PATH.")
    generate_umap_images(args.input, args.figures, args.catalog)


if __name__ == "__main__":
    main()
