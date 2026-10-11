UMAP_METADATA_COLUMNS = [
    "cell_complexity",
    "cell_type2",
    "celltypist_confidence",
    "celltypist_prediction",
    "dataset",
    "detected_genes",
    "diagnosis",
    "mito_perc",
    "n_genes_by_counts",
    "pct_counts_hb",
    "pct_counts_mt",
    "pct_counts_rp",
    "predicted_doublet",
]

UMAP_STATIC_IMAGE_FILES = {
    column: f"umap/{column}.png" for column in UMAP_METADATA_COLUMNS
}
