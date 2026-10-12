############################################################
# Cell-type annotation from canonical marker signatures
# Dataset: GSE155698 (PDAC)
#
# Clusters are labelled from data, not from a fixed position list:
# each canonical marker set is scored per cluster (AddModuleScore)
# and every cluster is assigned the cell type with the highest
# mean score. Clusters whose best score is not positive are left
# as "Unassigned".
############################################################

library(Seurat)

# Locate the repository root so the script runs from any working directory.
local({
  a <- commandArgs(trailingOnly = FALSE)
  f <- sub("^--file=", "", a[grep("^--file=", a)])
  root <- if (length(f)) dirname(dirname(normalizePath(f)))
          else if (file.exists("scripts")) normalizePath(".")
          else normalizePath("..")
  setwd(root)
})

set.seed(123)

#-----------------------------
# Load clustered object
#-----------------------------
pdac <- readRDS("data/processed/PDAC_UMAP.rds")

# Collapse the per-sample count layers into one (Seurat v5).
pdac <- JoinLayers(pdac)

#-----------------------------
# Canonical marker signatures
#-----------------------------
signatures <- list(

  Acinar = c("PRSS1", "CPA1", "CPA2", "CELA3A", "CELA3B", "CTRB1",
             "CTRB2", "CTRC", "CLPS", "PNLIP", "AMY2A", "AMY2B",
             "CEL", "REG1A", "REG1B", "SYCN", "GP2"),

  Ductal_Epithelial = c("KRT19", "KRT8", "KRT18", "EPCAM", "KRT7",
                        "MUC1", "CEACAM6", "TFF1", "TFF2", "S100P",
                        "LGALS4", "KRT17", "SOX9"),

  Endocrine = c("CHGA", "CHGB", "INS", "GCG", "SST", "NEUROD1",
                "PCSK1", "PCSK2"),

  T_Cell = c("CD3D", "CD3E", "CD3G", "CD2", "IL7R", "TRAC",
             "CD247", "LTB"),

  CD8_T = c("CD8A", "CD8B", "GZMK", "CCL5", "GZMA"),

  NK = c("NKG7", "GNLY", "KLRD1", "NCAM1", "PRF1", "KLRF1",
         "KLRB1"),

  B_Cell = c("CD79A", "CD79B", "MS4A1", "CD19", "IGKC", "IGHM",
             "BANK1"),

  Plasma_Cell = c("MZB1", "JCHAIN", "IGHG1", "IGHA1", "DERL3",
                  "XBP1", "TNFRSF17"),

  Monocyte = c("FCN1", "VCAN", "CD14", "LYZ", "S100A12", "CLEC7A",
               "MNDA"),

  Macrophage = c("CD68", "CD163", "MRC1", "MSR1", "C1QA", "C1QB",
                 "APOE", "FOLR2"),

  Neutrophil = c("FCGR3B", "CSF3R", "CXCR2", "BASP1", "S100A8",
                 "S100A9", "IL1B"),

  Mast_Cell = c("TPSAB1", "TPSB2", "CPA3", "MS4A2", "KIT", "HDC"),

  Fibroblast_CAF = c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM",
                     "FAP", "PDGFRA", "POSTN"),

  Endothelial = c("PECAM1", "VWF", "CDH5", "CLDN5", "RAMP2",
                  "EGFL7"),

  Pericyte = c("RGS5", "NOTCH3", "PDGFRB", "ACTA2", "TAGLN",
               "MCAM"),

  Cycling = c("MKI67", "TOP2A", "CENPF", "STMN1", "PCNA", "CENPE",
              "UBE2C"),

  Hypoxic_Tumor = c("VEGFA", "NDRG1", "SLC2A1", "CA9", "LDHA",
                    "ADM", "HILPDA")

)

# Keep only genes present in the object.
present <- rownames(pdac)
n_missing <- sum(!unlist(signatures) %in% present)
signatures <- lapply(signatures, function(g) intersect(g, present))
cat("Signature genes absent from object:", n_missing, "\n")

#-----------------------------
# Score signatures per cell
#-----------------------------
pdac <- AddModuleScore(
  object = pdac,
  features = signatures,
  name = "sig",
  seed = 123
)

sig_cols <- paste0("sig", seq_along(signatures))

#-----------------------------
# Mean score per cluster
#-----------------------------
clusters <- levels(Idents(pdac))

score_mat <- do.call(rbind, lapply(clusters, function(cl) {
  cells <- WhichCells(pdac, idents = cl)
  colMeans(pdac@meta.data[cells, sig_cols, drop = FALSE])
}))
rownames(score_mat) <- clusters
colnames(score_mat) <- names(signatures)

best <- max.col(score_mat, ties.method = "first")
labels <- names(signatures)[best]
labels[score_mat[cbind(seq_along(clusters), best)] <= 0] <- "Unassigned"
names(labels) <- clusters

cat("\nCluster assignments:\n")
print(data.frame(cluster = clusters, celltype = labels,
                 best_score = round(score_mat[cbind(seq_along(clusters), best)], 4),
                 row.names = NULL))

#-----------------------------
# Apply labels
#-----------------------------
pdac <- RenameIdents(pdac, labels)
pdac$celltype <- Idents(pdac)

#-----------------------------
# Save score evidence
#-----------------------------
dir.create("results/Tables", showWarnings = FALSE, recursive = TRUE)
write.csv(
  data.frame(cluster = clusters, score_mat, celltype = labels,
             check.names = FALSE, row.names = NULL),
  "results/Tables/CellType_Signature_Scores.csv",
  row.names = FALSE
)

#-----------------------------
# Save object
#-----------------------------
saveRDS(
  pdac,
  file = "data/processed/PDAC_Annotated.rds"
)

#-----------------------------
# Show cell counts
#-----------------------------
print(table(pdac$celltype))

cat("\n==============================\n")
cat("Cell annotation completed!\n")
cat("Annotated object saved to:\n")
cat("data/processed/PDAC_Annotated.rds\n")
cat("==============================\n")
