############################################################
# 08_DEG_Tumor_Normal.R
# Tumor vs Normal Differential Expression Analysis
# PDAC scRNA-seq
############################################################

library(Seurat)
library(dplyr)
library(ggplot2)

############################################################
# Load data
############################################################

pdac <- readRDS(
"data/processed/PDAC_Annotated.rds"
)

DefaultAssay(pdac) <- "RNA"

# Fix Seurat v5 multi-layer issue
pdac[["RNA"]] <- JoinLayers(pdac[["RNA"]])

Idents(pdac) <- "condition"

############################################################
# Output folders
############################################################

dir.create(
"results/DEG",
recursive = TRUE,
showWarnings = FALSE
)


############################################################
# Differential expression
############################################################

deg <- FindMarkers(
pdac,
ident.1 = "Tumor",
ident.2 = "Normal",
test.use = "wilcox",
min.pct = 0.25,
logfc.threshold = 0.25
)


############################################################
# Save DEG table
############################################################

deg$gene <- rownames(deg)

deg <- deg %>%
arrange(p_val_adj)


write.csv(
deg,
"results/DEG/Tumor_vs_Normal_DEG.csv",
row.names = FALSE
)


############################################################
# Significant DEGs
############################################################

sig_deg <- deg %>%
filter(
p_val_adj < 0.05,
abs(avg_log2FC) > 0.5
)


write.csv(
sig_deg,
"results/DEG/Significant_Tumor_vs_Normal_DEG.csv",
row.names = FALSE
)


############################################################
# Volcano plot
############################################################


volcano <- ggplot(
deg,
aes(
x = avg_log2FC,
y = -log10(p_val_adj)
)
) +

geom_point(
alpha = 0.5,
size = 1
) +

geom_point(
data = sig_deg,
aes(
x = avg_log2FC,
y = -log10(p_val_adj)
),
size = 1.5
) +

theme_classic(base_size = 15) +

labs(
title="Tumor vs Normal Differential Expression",
x="Log2 Fold Change",
y="-Log10 Adjusted P-value"
)


ggsave(
"results/DEG/Tumor_vs_Normal_Volcano.png",
volcano,
width=8,
height=6,
dpi=600,
bg="white"
)


############################################################
# Top DEG heatmap
############################################################

top_genes <- sig_deg %>%
arrange(desc(abs(avg_log2FC))) %>%
head(30) %>%
pull(gene)


pdac <- ScaleData(
pdac,
features = top_genes,
verbose = FALSE
)


heat <- DoHeatmap(
pdac,
features = top_genes,
group.by="condition",
size=4
) +
ggtitle(
"Top Differentially Expressed Genes: Tumor vs Normal"
)


ggsave(
"results/DEG/Tumor_Normal_DEG_Heatmap.png",
heat,
width=10,
height=10,
dpi=600,
bg="white"
)


############################################################

cat("
=====================================
Tumor vs Normal DEG Completed

Results:
results/DEG/

=====================================
")
