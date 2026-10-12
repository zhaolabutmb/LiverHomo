Static UMAP generation
======================

From the project root, run::

    .venv/bin/python scripts/generate_umap_images.py

The generator reads ``DATA_PATH`` from the environment or the shared ``.env``
(see ``.env.example``; override its location with ``ENV_FILE``).
Use ``--input`` to specify another AnnData file. It reads only the saved
``obsm/X_umap`` coordinates and observation columns, not expression matrices.

One PNG is written to ``figs/umap/`` for each ``UMAP_METADATA_COLUMNS`` entry
in ``umap_config.py``. Every image includes all cells, cell types and diseases,
without subsampling. Numeric metadata uses a colorbar; categorical metadata
uses a legend. Missing color values are retained, not filtered out.

``figs/umap/manifest.json`` records cell counts and coloring information.
After all images are generated, the generator updates only ``umap_columns``
in the pregenerated ``static/filter_catalog.json``.

The Atlas metadata dropdown changes coloring only. There are no UMAP disease
or cell-type filters. Restart Gunicorn workers after changes to the shared
UMAP configuration; regenerating images alone does not require a restart.
