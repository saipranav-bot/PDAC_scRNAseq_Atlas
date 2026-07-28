############################################################
# DEG_Top15_DotPlot.R
# Publication Quality Top 15 DEG DotPlot
############################################################

library(Seurat)
library(ggplot2)
library(dplyr)

############################################################
# Load data
############################################################

pdac <- readRDS("data/processed/PDAC_Annotated.rds")

DefaultAssay(pdac) <- "RNA"

# Seurat v5 fix
pdac[["RNA"]] <- JoinLayers(pdac[["RNA"]])

############################################################
# Load DEG results
############################################################

deg <- read.csv(
  "results/DEG/Tumor_vs_Normal_DEG.csv"
)

############################################################
# Select significant genes
############################################################

deg_clean <- deg %>%
  filter(
    p_val_adj < 0.05,
    abs(avg_log2FC) > 0.5
  )

# Remove housekeeping/background genes

deg_clean <- deg_clean %>%
  filter(
    !grepl(
      "^MT-|^RPL|^RPS|^LINC|^RP11",
      gene,
      ignore.case = TRUE
    )
  )


############################################################
# Select balanced Top 15 genes
############################################################

top_tumor <- deg_clean %>%
  filter(avg_log2FC > 0) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 8)


top_normal <- deg_clean %>%
  filter(avg_log2FC < 0) %>%
  arrange(p_val_adj) %>%
  slice_head(n = 7)


top15 <- bind_rows(
  top_tumor,
  top_normal
)


genes <- top15$gene


write.csv(
  top15,
  "results/DEG/Top15_DEG_Selected.csv",
  row.names = FALSE
)


############################################################
# Generate DotPlot
############################################################

dot <- DotPlot(
  pdac,
  features = genes,
  group.by = "condition",
  dot.scale = 7
) +
  
  RotatedAxis() +
  
  scale_color_gradient(
    low = "lightgrey",
    high = "red"
  ) +
  
  ggtitle(
    "Top 15 Differentially Expressed Genes\nTumor vs Normal"
  ) +
  
  theme_classic(base_size = 15) +
  
  theme(
    
    plot.title = element_text(
      hjust = 0.5,
      face = "bold",
      size = 16
    ),
    
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 11,
      face = "italic"
    ),
    
    axis.text.y = element_text(
      size = 13,
      face = "bold"
    ),
    
    axis.title = element_blank(),
    
    legend.position = "right"
  )


############################################################
# Save
############################################################


ggsave(
  "results/DEG/Top15_DEG_DotPlot.png",
  dot,
  width = 9,
  height = 6,
  dpi = 600,
  bg = "white"
)


cat(
"
=========================================
Top15 DEG DotPlot completed

Output:
results/DEG/Top15_DEG_DotPlot.png
results/DEG/Top15_DEG_Selected.csv

=========================================
"
)
