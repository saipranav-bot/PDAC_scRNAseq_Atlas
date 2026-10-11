# PDAC scRNA-seq Atlas

Single-cell RNA-seq analysis of pancreatic ductal adenocarcinoma (PDAC),
built with **Seurat** (R). Cells from tumour and adjacent normal pancreas are
integrated, quality-controlled, clustered, annotated by canonical markers, and
compared (tumour vs normal) with differential expression, GO/KEGG enrichment and
GSEA.

## Pipeline

Run the numbered scripts in order from the **repository root**:

| Step | Script | Purpose |
|------|--------|---------|
| 1 | `scripts/01_load_and_merge_samples.R` | Load and merge samples |
| 2 | `scripts/02_QC_analysis.R` | Quality control and filtering |
| 3 | `scripts/03_Normalization_PCA.R` | Normalisation, HVGs, PCA |
| 4 | `scripts/04_Clustering_UMAP.R` | Clustering and UMAP embedding |
| 5 | `scripts/05_Find_Markers.R` | Cluster marker genes |
| 6 | `scripts/06_Cell_Annotation.R` | Cell-type annotation |
| 7 | `scripts/07_Visualization.R` | Publication visualisations |
| 8 | `scripts/08_DEG_Tumor_Normal.R` | Tumour vs normal DEGs |
| 9 | `scripts/09_Pathway_Enrichment.R` | GO and KEGG enrichment |
| 10 | `scripts/14_GSEA_Analysis.R` | Gene-set enrichment analysis |

Publication figures are produced by `scripts/DEG_Publication_Final.R`,
`scripts/DEG_Publication_Plots.R`, `scripts/DEG_Top15_DotPlot.R` and
`scripts/Top15_DEG_Plot.R`. `scripts/export_markers.py` and
`scripts/filter_markers.py` are scanpy helpers for exporting/filtering marker
tables.

Each script locates the repository root automatically, so it can be run from any
working directory.

## Repository structure

```text
PDAC_scRNAseq_Atlas/
├── scripts/       # analysis scripts (R + two Python helpers)
├── data/
│   ├── raw/       # raw 10x matrices (not tracked)
│   └── processed/ # Seurat/scanpy objects (not tracked)
├── results/       # tables, markers, pathways, clustering, QC
├── figures/       # summary figures
└── README.md
```

Raw and processed matrices are not committed; fetch the source data from GEO
before running the pipeline:

```bash
bash scripts/00_download_data.sh
```

This downloads the ten GSE155698 samples used here (7 tumour + 3 adjacent-normal)
into `data/raw/<sample>/filtered_feature_bc_matrix/`. Note that `PDAC_TISSUE_16`
is deposited by GEO in the older Cell Ranger `filtered_gene_bc_matrices` layout;
`01_load_and_merge_samples.R` handles both formats automatically.
