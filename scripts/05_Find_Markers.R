############################################################
# Find Marker Genes for Each Cluster (Seurat v5)
############################################################

library(Seurat)
library(dplyr)

setwd("~/Cancer_scRNA_Project")

############################################################
# Load clustered object
############################################################

pdac <- readRDS(
  "data/processed/PDAC_UMAP.rds"
)

############################################################
# Join layers (Required for Seurat v5)
############################################################

pdac <- JoinLayers(pdac)

############################################################
# Check clusters
############################################################

print(table(Idents(pdac)))

############################################################
# Find marker genes
############################################################

markers <- FindAllMarkers(
  object = pdac,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25,
  verbose = TRUE
)

############################################################
# Check results
############################################################

cat("\nNumber of marker genes found:\n")
print(nrow(markers))

print(head(markers))

############################################################
# Save complete marker table
############################################################

write.csv(
  markers,
  "results/Cluster_Markers.csv",
  row.names = FALSE
)

############################################################
# Top 10 markers per cluster
############################################################

top10 <- markers %>%
  group_by(cluster) %>%
  slice_max(
    order_by = avg_log2FC,
    n = 10,
    with_ties = FALSE
  )

############################################################
# Save top markers
############################################################

write.csv(
  top10,
  "results/Top10_Markers_Per_Cluster.csv",
  row.names = FALSE
)

############################################################
# Print top markers
############################################################

print(top10)

############################################################
# Save object
############################################################

saveRDS(
  pdac,
  "data/processed/PDAC_Markers.rds"
)

cat("\n=====================================\n")
cat("Marker analysis completed successfully\n")
cat("Files created:\n")
cat("results/Cluster_Markers.csv\n")
cat("results/Top10_Markers_Per_Cluster.csv\n")
cat("=====================================\n")
