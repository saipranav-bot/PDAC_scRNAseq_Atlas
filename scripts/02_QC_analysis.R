############################################################
# Single Cell RNA-seq QC Analysis
# Dataset: GSE155698 PDAC
############################################################


##############################
# Load libraries
##############################

library(Seurat)
library(dplyr)
library(ggplot2)


##############################
# Set directory
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


##############################
# Load merged object
##############################

pdac_combined <- readRDS(
  "data/processed/PDAC_merged_raw.rds"
)


##############################
# Calculate mitochondrial genes
##############################

pdac_combined[["percent.mt"]] <- PercentageFeatureSet(
  pdac_combined,
  pattern = "^MT-"
)


##############################
# View metadata
##############################

head(
  pdac_combined@meta.data
)


##############################
# QC summary
##############################

summary(
  pdac_combined$nFeature_RNA
)

summary(
  pdac_combined$nCount_RNA
)

summary(
  pdac_combined$percent.mt
)



##############################
# QC violin plots
##############################

pdf(
  "results/QC_violin_before_filtering.pdf",
  width = 10,
  height = 5
)


VlnPlot(
  pdac_combined,
  features =
    c(
      "nFeature_RNA",
      "nCount_RNA",
      "percent.mt"
    ),
  ncol = 3,
  pt.size = 0
)


dev.off()



##############################
# QC scatter plots
##############################

pdf(
  "results/QC_scatter_before_filtering.pdf",
  width = 8,
  height = 4
)


FeatureScatter(
  pdac_combined,
  feature1 = "nCount_RNA",
  feature2 = "nFeature_RNA"
)


FeatureScatter(
  pdac_combined,
  feature1 = "nCount_RNA",
  feature2 = "percent.mt"
)


dev.off()



##############################
# Filtering cells
##############################

pdac_filtered <- subset(

  pdac_combined,

  subset =
    nFeature_RNA > 200 &
    nFeature_RNA < 6000 &
    percent.mt < 20

)



##############################
# Compare cell numbers
##############################

print(
  paste(
    "Before filtering:",
    ncol(pdac_combined)
  )
)


print(
  paste(
    "After filtering:",
    ncol(pdac_filtered)
  )
)



##############################
# Save QC filtered object
##############################

saveRDS(

  pdac_filtered,

  file =
    "data/processed/PDAC_filtered_QC.rds"

)



print(
  "QC filtering completed"
)
