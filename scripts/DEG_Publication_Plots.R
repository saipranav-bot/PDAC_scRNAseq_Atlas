############################################################
# DEG_Publication_Plots.R
# Publication Quality Tumor vs Normal DEG Visualization
############################################################


library(Seurat)
library(SeuratObject)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(pheatmap)


set.seed(123)


dir.create(
  "results/DEG",
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# Load Seurat object
############################################################

pdac <- readRDS(
  "data/processed/PDAC_Annotated.rds"
)


DefaultAssay(pdac) <- "RNA"


# Seurat v5 layer correction
pdac[["RNA"]] <- JoinLayers(
  pdac[["RNA"]]
)



############################################################
# Load DEG results
############################################################

deg <- read.csv(
  "results/DEG/Tumor_vs_Normal_DEG.csv"
)



############################################################
# Volcano Plot
############################################################


deg <- deg %>%
  mutate(
    Significance = case_when(

      p_val_adj < 0.05 &
        avg_log2FC > 0.5 ~ "Tumor",

      p_val_adj < 0.05 &
        avg_log2FC < -0.5 ~ "Normal",

      TRUE ~ "NS"
    )
  )



# Select top 15 most significant genes

label_genes <- deg %>%
  filter(Significance != "NS") %>%
  arrange(p_val_adj) %>%
  slice_head(n = 15)



volcano <- ggplot(
  deg,
  aes(
    x = avg_log2FC,
    y = -log10(p_val_adj),
    color = Significance
  )
) +

geom_point(
  size = 2.2,
  alpha = 0.75
) +

scale_color_manual(
  values=c(
    "Tumor"="red",
    "Normal"="blue",
    "NS"="grey75"
  )
) +


geom_text_repel(
  data = label_genes,
  aes(label=gene),
  size=4,
  max.overlaps=15,
  box.padding=0.5
) +


theme_classic(
  base_size = 16
) +

theme(
  plot.title =
    element_text(
      hjust=0.5,
      face="bold"
    ),

  legend.position="right"
) +


labs(
  title="Differential Gene Expression\nTumor vs Normal PDAC",
  x="log2 Fold Change",
  y="-log10 Adjusted P-value"
)



ggsave(
  "results/DEG/Tumor_Normal_Volcano_Publication.png",
  volcano,
  width=10,
  height=8,
  dpi=600,
  bg="white"
)




############################################################
# Top 15 DEG Heatmap
############################################################


# Select top 15 strongest DEGs

top15 <- deg %>%
  filter(
    p_val_adj < 0.05
  ) %>%
  arrange(
    p_val_adj,
    desc(abs(avg_log2FC))
  ) %>%
  slice_head(n=15)



write.csv(
  top15,
  "results/DEG/Top15_DEG_Selected.csv",
  row.names=FALSE
)



genes <- top15$gene



############################################################
# Average expression matrix
############################################################


avg_expression <- AverageExpression(
  pdac,
  assays="RNA",
  features=genes,
  group.by="condition"
)$RNA



# scale genes

heat_matrix <- t(
  scale(
    t(avg_expression)
  )
)



############################################################
# Column annotation
############################################################


annotation <- data.frame(
  Condition = colnames(heat_matrix)
)


rownames(annotation) <- colnames(heat_matrix)



############################################################
# Save heatmap
############################################################


png(
  "results/DEG/Tumor_Normal_Top15_DEG_Heatmap.png",
  width=2800,
  height=2200,
  res=300
)


pheatmap(
  heat_matrix,

  annotation_col =
    annotation,

  cluster_rows=TRUE,

  cluster_cols=FALSE,

  color =
    colorRampPalette(
      c(
        "blue",
        "white",
        "red"
      )
    )(100),

  fontsize_row=14,

  fontsize_col=16,

  angle_col=0,

  border_color=NA,

  main=
    "Top 15 Differentially Expressed Genes\nTumor vs Normal"
)


dev.off()



cat(
"
=========================================
Publication DEG plots completed
=========================================

Generated:

1. Tumor_Normal_Volcano_Publication.png

2. Tumor_Normal_Top15_DEG_Heatmap.png

3. Top15_DEG_Selected.csv

=========================================
"
)
