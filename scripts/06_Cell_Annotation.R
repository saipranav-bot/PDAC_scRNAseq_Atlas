library(Seurat)

#-----------------------------
# Load clustered object
#-----------------------------
pdac <- readRDS("data/processed/PDAC_UMAP.rds")
#-----------------------------
# Cell type annotation
#-----------------------------

celltypes <- c(

"Acinar",
"CD8_T",
"Tumor_Epithelial",
"Naive_T",
"Monocyte",
"Neutrophil",
"TAM",
"Resident_Macrophage",
"NK",
"Inflammatory_Monocyte",
"Tumor_Epithelial_2",
"Mucinous_Tumor",
"Hypoxic_Tumor",
"Exhausted_CD8",
"Gastric_like_Tumor",
"B_Cell",
"Pericyte",
"CAF",
"T_Cell",
"Acinar_like",
"Endothelial",
"Mast_Cell",
"Rare_Epithelial",
"Cycling_Cells",
"Endocrine",
"Tumor_Epithelial_3",
"Cycling_Tumor"

)

names(celltypes) <- levels(pdac)

pdac <- RenameIdents(pdac, celltypes)

pdac$celltype <- Idents(pdac)

#-----------------------------
# Save object
#-----------------------------

saveRDS(
    pdac,
    file="data/processed/PDAC_Annotated.rds"
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
