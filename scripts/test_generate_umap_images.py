import json
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import anndata as ad
import numpy as np
import pandas as pd
from PIL import Image

from scripts.generate_umap_images import generate_umap_images
from umap_config import UMAP_METADATA_COLUMNS, UMAP_STATIC_IMAGE_FILES


class GenerateUmapImagesTests(unittest.TestCase):
    def setUp(self):
        self.directory = tempfile.TemporaryDirectory()
        self.addCleanup(self.directory.cleanup)
        self.root = Path(self.directory.name)
        self.input_path = self.root / "source.h5ad"
        self.catalog_path = self.root / "filter_catalog.json"
        self.catalog_path.write_text(json.dumps({"umap_columns": [], "metadata": {"keep": True}}))
        observations = pd.DataFrame({
            column: np.array([0.0, 1.0, 2.0, np.nan])
            for column in UMAP_METADATA_COLUMNS
        }, index=["c1", "c2", "c3", "c4"])
        for column in ("cell_type2", "diagnosis", "dataset", "celltypist_prediction"):
            observations[column] = pd.Categorical(["A", "B", "A", None])
        observations["predicted_doublet"] = pd.Categorical([False, True, False, True])
        self.adata = ad.AnnData(obs=observations)
        self.adata.obsm["X_umap"] = np.array([[0, 0], [1, 1], [2, 0], [3, 1]], dtype=float)

    def test_all_columns_and_cells_are_generated(self):
        self.adata.write_h5ad(self.input_path)
        manifest = generate_umap_images(self.input_path, self.root / "figs", self.catalog_path)
        self.assertFalse(manifest["subsampled"])
        self.assertEqual(manifest["filters"], {})
        self.assertEqual(manifest["cells"], 4)
        self.assertEqual([plot["column"] for plot in manifest["plots"]], UMAP_METADATA_COLUMNS)
        for plot in manifest["plots"]:
            self.assertEqual(plot["cells_plotted"], 4)
            with Image.open(self.root / "figs" / plot["file"]) as image:
                self.assertEqual(image.format, "PNG")
                image.verify()
        self.assertEqual(
            json.loads(self.catalog_path.read_text()),
            {"umap_columns": UMAP_METADATA_COLUMNS, "metadata": {"keep": True}},
        )
        numeric = manifest["plots"][0]
        self.assertEqual(numeric["color_type"], "numeric")
        self.assertEqual(numeric["missing_values"], 1)
        categorical = next(plot for plot in manifest["plots"] if plot["column"] == "diagnosis")
        self.assertIn("Unknown", categorical["categories"])

    def test_missing_column_fails_without_updating_catalog(self):
        del self.adata.obs["diagnosis"]
        self.adata.write_h5ad(self.input_path)
        with self.assertRaisesRegex(ValueError, "diagnosis"):
            generate_umap_images(self.input_path, self.root / "figs", self.catalog_path)
        self.assertEqual(json.loads(self.catalog_path.read_text())["umap_columns"], [])
        self.assertFalse((self.root / "figs").exists())

    def test_invalid_coordinates_fail_without_updating_catalog(self):
        self.adata.obsm["X_umap"][0, 0] = np.nan
        self.adata.write_h5ad(self.input_path)
        with self.assertRaisesRegex(ValueError, "finite"):
            generate_umap_images(self.input_path, self.root / "figs", self.catalog_path)
        self.assertEqual(json.loads(self.catalog_path.read_text())["umap_columns"], [])

    def test_render_failure_does_not_publish_catalog(self):
        self.adata.write_h5ad(self.input_path)
        with patch("scripts.generate_umap_images.render_umap", side_effect=OSError("disk failure")):
            with self.assertRaisesRegex(OSError, "disk failure"):
                generate_umap_images(self.input_path, self.root / "figs", self.catalog_path)
        self.assertEqual(json.loads(self.catalog_path.read_text())["umap_columns"], [])

    def test_static_umap_routes_and_metadata_selector(self):
        import app as app_module

        self.adata.write_h5ad(self.input_path)
        generate_umap_images(self.input_path, self.root / "figs", self.catalog_path)
        with patch.object(app_module.app, "root_path", str(self.root)), patch.object(
            app_module, "FILTER_CATALOG_PATH", str(self.catalog_path)
        ):
            client = app_module.app.test_client()
            for column in UMAP_METADATA_COLUMNS:
                response = client.get(f"/umap_image/{column}.png")
                self.assertEqual(response.status_code, 200)
                self.assertEqual(response.mimetype, "image/png")
                response.close()
                metadata_response = client.get(f"/api/umap_metadata?column={column}")
                self.assertEqual(metadata_response.status_code, 200)
                self.assertEqual(metadata_response.json["column"], column)
            self.assertEqual(client.get("/umap_image/unknown.png").status_code, 404)
            self.assertEqual(client.get("/api/umap_metadata?column=unknown").status_code, 404)
        with patch.object(app_module, "get_static_umap_metadata_columns", return_value=UMAP_METADATA_COLUMNS), patch.object(
            app_module, "load_filter_catalog", return_value={"metadata": {}}
        ):
            html = client.get("/atlas").get_data(as_text=True)
            self.assertIn('/umap_image/diagnosis.png', html)
            self.assertIn('id="umap-column"', html)
            for column in UMAP_METADATA_COLUMNS:
                self.assertIn(f'value="{column}"', html)
            self.assertIn("Select a metadata column", html)
        self.assertEqual(set(UMAP_STATIC_IMAGE_FILES), set(UMAP_METADATA_COLUMNS))


if __name__ == "__main__":
    unittest.main()
