############################################################
# Clustering and UMAP
############################################################

library(Seurat)
library(ggplot2)

setwd("~/Cancer_scRNA_Project")

##############################
# Load PCA object
##############################

pdac <- readRDS(
  "data/processed/PDAC_PCA.rds"
)

##############################
# Find neighbors
##############################

pdac <- FindNeighbors(
  pdac,
  dims = 1:20
)

##############################
# Cluster cells
##############################

pdac <- FindClusters(
  pdac,
  resolution = 0.5
)

##############################
# Run UMAP
##############################

pdac <- RunUMAP(
  pdac,
  dims = 1:20
)

##############################
# Show cluster sizes
##############################

print(table(pdac$seurat_clusters))

##############################
# UMAP by cluster
##############################

pdf(
  "results/UMAP_clusters.pdf",
  width = 8,
  height = 6
)

print(
  DimPlot(
    pdac,
    reduction = "umap",
    label = TRUE
  )
)

dev.off()

##############################
# UMAP by condition
##############################

pdf(
  "results/UMAP_condition.pdf",
  width = 8,
  height = 6
)

print(
  DimPlot(
    pdac,
    reduction = "umap",
    group.by = "condition"
  )
)

dev.off()

##############################
# UMAP by sample
##############################

pdf(
  "results/UMAP_sample.pdf",
  width = 10,
  height = 8
)

print(
  DimPlot(
    pdac,
    reduction = "umap",
    group.by = "sample"
  )
)

dev.off()

##############################
# Save object
##############################

saveRDS(
  pdac,
  "data/processed/PDAC_UMAP.rds"
)

cat("\nUMAP analysis completed successfully.\n")
