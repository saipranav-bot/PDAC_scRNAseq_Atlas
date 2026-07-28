############################################################
# Top15_DEG_Plot.R
# Publication Quality Top 15 Tumor vs Normal DEG Visualization
############################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(pheatmap)

############################################################
# Create output folder
############################################################

dir.create(
  "results/DEG/Top15",
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# Load DEG results
############################################################

deg <- read.csv(
  "results/DEG/Significant_Tumor_vs_Normal_DEG.csv"
)


############################################################
# Select Top 15 Significant DEGs
############################################################

top15 <- deg %>%
  filter(!is.na(p_val_adj)) %>%
  arrange(
    p_val_adj,
    desc(abs(avg_log2FC))
  ) %>%
  slice_head(n = 15)


write.csv(
  top15,
  "results/DEG/Top15/Top15_DEG_Genes.csv",
  row.names = FALSE
)


############################################################
# Volcano Plot with Top 15 Labels
############################################################

deg <- deg %>%
  mutate(
    Significance = case_when(
      p_val_adj < 0.05 &
        abs(avg_log2FC) > 1 ~ "Significant",
      TRUE ~ "Not Significant"
    )
  )


top15 <- top15 %>%
  mutate(
    Significance = "Significant"
  )


volcano <- ggplot(
  deg,
  aes(
    x = avg_log2FC,
    y = -log10(p_val_adj),
    color = Significance
  )
) +

geom_point(
  size = 1.8,
  alpha = 0.6
) +

geom_text_repel(
  data = top15,
  aes(
    label = gene
  ),
  size = 4,
  max.overlaps = 20,
  box.padding = 0.5
) +

theme_classic(
  base_size = 15
) +

labs(
  title = "Tumor vs Normal Differential Expression",
  subtitle = "Top 15 Significant Genes",
  x = "log2 Fold Change",
  y = "-log10 Adjusted P-value"
) +

theme(
  plot.title = element_text(
    hjust = 0.5,
    face = "bold",
    size = 18
  ),
  plot.subtitle = element_text(
    hjust = 0.5
  ),
  legend.position = "top"
)


ggsave(
  "results/DEG/Top15/Tumor_Normal_Volcano_Top15.png",
  volcano,
  width = 10,
  height = 8,
  dpi = 600,
  bg = "white"
)


############################################################
# Heatmap Top 15 DEG
############################################################

pdac <- readRDS(
  "data/processed/PDAC_Annotated.rds"
)

DefaultAssay(pdac) <- "RNA"


############################################################
# Join Seurat v5 layers
############################################################

pdac[["RNA"]] <- JoinLayers(
  pdac[["RNA"]]
)


############################################################
# Extract genes present in object
############################################################

genes <- top15$gene

genes <- genes[
  genes %in% rownames(pdac)
]


############################################################
# Average expression Tumor vs Normal
############################################################

avg_exp <- AverageExpression(
  pdac,
  assays = "RNA",
  features = genes,
  group.by = "condition"
)$RNA


avg_exp <- as.matrix(avg_exp)


############################################################
# Save Heatmap
############################################################

png(
  "results/DEG/Top15/Tumor_Normal_Heatmap_Top15.png",
  width = 2600,
  height = 1800,
  res = 300
)


pheatmap(
  avg_exp,
  scale = "row",
  cluster_rows = TRUE,
  cluster_cols = FALSE,
  border_color = NA,
  fontsize_row = 12,
  fontsize_col = 14,
  main = "Top 15 Differentially Expressed Genes"
)


dev.off()


############################################################
# Completion message
############################################################

cat("
=====================================
Top 15 DEG Visualization Completed

Files saved:
results/DEG/Top15/

- Top15_DEG_Genes.csv
- Tumor_Normal_Volcano_Top15.png
- Tumor_Normal_Heatmap_Top15.png

=====================================
")
