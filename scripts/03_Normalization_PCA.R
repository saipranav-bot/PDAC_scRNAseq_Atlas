############################################################
# Normalization, Variable Features, Scaling and PCA
# Dataset: GSE155698 (PDAC)
############################################################

##############################
# Load libraries
##############################

library(Seurat)
library(ggplot2)

##############################
# Set working directory
##############################

# Locate the repository root so the script runs from any working directory.
local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  root <- if (length(f)) dirname(dirname(normalizePath(f)))
          else if (file.exists("scripts")) normalizePath(".")
          else normalizePath("..")
  setwd(root)
})

# Fixed seed so PCA (randomised SVD) is reproducible across runs.
set.seed(123)

##############################
# Load QC-filtered object
##############################

pdac <- readRDS(
  "data/processed/PDAC_filtered_QC.rds"
)

##############################
# Normalize data
##############################

pdac <- NormalizeData(
  pdac,
  normalization.method = "LogNormalize",
  scale.factor = 10000,
  verbose = TRUE
)

##############################
# Find highly variable genes
##############################

pdac <- FindVariableFeatures(
  pdac,
  selection.method = "vst",
  nfeatures = 2000,
  verbose = TRUE
)

##############################
# Top variable genes
##############################

top20 <- head(
  VariableFeatures(pdac),
  20
)

cat("\nTop 20 Variable Genes:\n")
print(top20)

##############################
# Plot variable features
##############################

pdf(
  "results/Variable_Features.pdf",
  width = 8,
  height = 6
)

plot1 <- VariableFeaturePlot(pdac)

plot2 <- LabelPoints(
  plot = plot1,
  points = top20,
  repel = TRUE
)

print(plot2)

dev.off()

##############################
# Scale ONLY variable genes
##############################

pdac <- ScaleData(
  pdac,
  features = VariableFeatures(pdac),
  verbose = TRUE
)

##############################
# Run PCA
##############################

pdac <- RunPCA(
  pdac,
  features = VariableFeatures(pdac),
  verbose = TRUE
)

##############################
# Print PCA summary
##############################

print(pdac[["pca"]])

##############################
# PCA plot
##############################

pdf(
  "results/PCA.pdf",
  width = 8,
  height = 6
)

print(
  DimPlot(
    pdac,
    reduction = "pca"
  )
)

dev.off()

##############################
# Elbow plot
##############################

pdf(
  "results/ElbowPlot.pdf",
  width = 6,
  height = 5
)

print(
  ElbowPlot(
    pdac,
    ndims = 50
  )
)

dev.off()

##############################
# Save processed object
##############################

saveRDS(
  pdac,
  file = "data/processed/PDAC_PCA.rds"
)

cat("\n=====================================\n")
cat("Normalization completed successfully\n")
cat("PCA completed successfully\n")
cat("Output saved:\n")
cat("data/processed/PDAC_PCA.rds\n")
cat("=====================================\n")
