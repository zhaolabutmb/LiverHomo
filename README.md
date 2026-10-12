# LiverHomo

**A single-cell atlas of the human liver in health and disease.**

LiverHomo integrates 16 public scRNA-seq datasets (2019–2024) into a harmonized atlas of
**1.16 million cells** across **26 cell types** and **12 groups**: Healthy, Adjacent Normal, Fetal,
and nine liver diseases. It connects gene expression and disease-associated changes with
cell–cell communication and context-specific protein interaction networks.

🔗 **Explore the atlas: [www.liverhomo.org](https://www.liverhomo.org)**

This repository contains the source code of the web portal.

## What you can do

| View | What it shows |
|---|---|
| **Atlas Explorer** | UMAPs of all cells coloured by metadata (cell type, disease, dataset, QC metrics) |
| **Gene Expression** | Expression of any of 31,351 genes by cell type and disease (histograms, violin plots) |
| **DEGs** | Differentially expressed genes for each disease and cell type versus Healthy, with volcano plots and downloadable tables |
| **Protein Interaction Network** | Context-specific protein–protein interaction (PPI) networks per disease and cell type |
| **EGO:BP Network** | Gene Ontology biological-process enrichment of the networks |
| **Cell–cell Interaction** | Interaction counts between cell types for each disease |
| **Ask LiverHomo** | An AI assistant that answers questions such as *"Is SPP1 upregulated in HCC macrophages?"* from the atlas data |

### Disease groups

| Short name | Group | Short name | Group |
|---|---|---|---|
| Healthy | Healthy liver (reference) | MASH | Metabolic dysfunction-associated steatohepatitis |
| Adjacent Normal | Tumour-adjacent normal tissue | MLN | Metastatic lymph nodes |
| Fetal | Fetal liver | PBC | Primary biliary cholangitis |
| ALD | Alcoholic liver disease | PSC | Primary sclerosing cholangitis |
| Cirrhosis | Cirrhosis | PVTT | Portal vein tumour thrombus |
| HCC | Hepatocellular carcinoma | iCCA | Intrahepatic cholangiocarcinoma |

## Methods in brief

- **Data:** human liver scRNA-seq from GEO, SRA and Tabula Sapiens; 10x Genomics-compatible
  platforms only. The full dataset list is on the
  [Resources page](https://www.liverhomo.org/resources), where you can also download the cell
  metadata (all 1.16 million cells, 150 MB, `.h5ad`).
- **Processing:** CellRanger v7.1.0 (GRCh38) → Scanpy QC, normalization and doublet removal →
  Harmony batch correction → PCA, UMAP, Leiden clustering.
- **Annotation:** CellTypist, validated with marker genes from the literature and CellMarker;
  T-cell subsets with starCAT.
- **Differential expression:** memento, each disease against Healthy. DEG tables list
  significant genes only.
- **Cell–cell interaction:** CellPhoneDB consensus networks per disease.
- **Protein interaction networks:** PrePPI, built from cell-type marker genes and from genes
  upregulated in each disease; annotated with UniProt and Gene Ontology.

## Citation

If you use LiverHomo in your research, please cite:

> Shi, Ziyu, et al. *LIVERDCP: A Disease-Cell-Protein Framework for Multi-Scale Modeling of
> Disease Biology.* bioRxiv, 10 Aug. 2026.
> [doi:10.64898/2026.08.07.743628](https://doi.org/10.64898/2026.08.07.743628)

---

## Running the portal

The sections below are for developers who want to run or modify the site.

### Requirements

- Python 3.11 and the packages in `requirements.txt`
- The atlas data files (about 38 GB). They are **not** included in this repository because of
  their size. [Contact us](https://www.liverhomo.org/contact) for access.
- An OpenAI API key, for the AI assistant only
- About 2 GB of free RAM. The app reads expression data gene by gene from disk rather than
  loading whole matrices.

### Data layout

The app finds its data through paths set in a `.env` file:

| Variable | Contents |
|---|---|
| `DATA_PATH` | Full processed AnnData (`adata_reprocessed_gzip.h5ad`); input for the scripts below |
| `METADATA_PATH` | Cell metadata only (`adata_metadata.h5ad`) |
| `GENE_EXPRESSION_DATA_PATH` | Expression matrices for the gene viewer (`adata_gene_expression.h5ad`) |
| `DEG_ROOT` | DEG tables: `DEG_<disease>.csv` and `<disease>/DEG_<disease>_<cell type>.csv` |
| `GRAPH_ROOT` | PPI edge lists: `<disease>_<cell type>_querylist1.xlsx` |
| `FIGS_ROOT` | Pre-rendered figures; UMAP images in `umap/` |
| `VOLCANO_ROOT`, `EGO_ROOT`, `CELL_INTERACTION_ROOT` | Volcano plots, GO networks, interaction heatmaps and counts |
| `GENE_VECTOR_CACHE_DIR` | Writable folder for the per-gene cache |

### Setup

```bash
git clone git@github.com:zhaolabutmb/LiverHomo.git
cd LiverHomo
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt

cp .env.example ~/liverhomo.env      # then edit the paths and add OPENAI_API_KEY
export ENV_FILE=~/liverhomo.env      # not needed on the production server
```

`.env.example` lists every variable. Keep the real `.env` outside the repository; it holds
secrets and is ignored by git. Without `ENV_FILE`, the app reads `/home/administrator/.env`.

### Generate derived files

From a new `DATA_PATH`, build the smaller files the app reads:

```bash
python scripts/generate_metadata_adata.py          # -> METADATA_PATH
python scripts/generate_gene_expression_adata.py   # -> GENE_EXPRESSION_DATA_PATH
python scripts/generate_umap_images.py             # -> FIGS_ROOT/umap/*.png
python scripts/generate_filter_catalog.py          # -> static/filter_catalog.json
```

### Start the server

```bash
gunicorn --config gunicorn.conf.py wsgi:app
```

It listens on `0.0.0.0:5000` with 1 worker and 4 threads. Keep one worker on small servers:
each worker holds its own copy of the cached data.

### AI assistant

The assistant (`/api/query`) uses OpenAI's `gpt-5.4-nano` with tools that query the atlas
directly: DEG search by gene, ranked DEGs, expression summaries, cell counts and figures.
It labels general background knowledge separately from atlas results, and declines questions
outside the atlas's scope, including medical advice. Set `ASSISTANT_MODEL` to use another model.

To check answer quality after changing the assistant, run the test questions in
`scripts/assistant_eval_cases.json`:

```bash
python scripts/eval_assistant.py --label my-change
```

A full run costs about US$0.02 in API usage.

### Repository layout

```
app.py                 Flask app: pages, data APIs, AI assistant
assistant_data.py      DEG index and disease/cell-type name matching for the assistant
expression_store.py    Low-memory readers for the expression and metadata files
env_config.py          Loads the shared .env (override with ENV_FILE)
umap_config.py         Metadata columns that have UMAP images
gunicorn.conf.py       Server settings
templates/, static/    Web pages, JavaScript, CSS, images
scripts/               Data preparation and assistant evaluation
```

## Contact

Zhao Lab, The University of Texas Medical Branch:
[zhaolabutmb.github.io](https://zhaolabutmb.github.io) ·
[contact form](https://www.liverhomo.org/contact)

## Acknowledgements

Supported by [The Welch Foundation](https://welch1.org/) and the
[UT System STARS Program](https://www.utsystem.edu/offices/health-affairs/initiatives-and-departments/facility-science-and-technology-acquisition-and-retention-stars-program).
