############################################################
# 07_Visualization.R
# Final Publication Quality PDAC scRNA-seq Visualization
############################################################

library(Seurat)
library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(1234)

############################################################
# Output directory
############################################################

dir.create(
  "results/Figures",
  recursive = TRUE,
  showWarnings = FALSE
)


############################################################
# Load Annotated Object
############################################################

pdac <- readRDS(
  "data/processed/PDAC_Annotated.rds"
)

DefaultAssay(pdac) <- "RNA"

# Seurat v5 fix
pdac[["RNA"]] <- JoinLayers(pdac[["RNA"]])

Idents(pdac) <- "celltype"


############################################################
# Global Theme
############################################################

theme_set(
  theme_classic(base_size = 15) +
    theme(
      plot.title = element_text(
        hjust = 0.5,
        face = "bold",
        size = 18
      ),
      axis.text = element_text(size = 11),
      axis.title = element_text(size = 13),
      legend.text = element_text(size = 10),
      panel.background = element_rect(fill="white"),
      plot.background = element_rect(fill="white")
    )
)


############################################################
# Figure 1
# UMAP Cell Types
############################################################

p1 <- DimPlot(
  pdac,
  reduction="umap",
  group.by="celltype",
  label=TRUE,
  repel=TRUE,
  pt.size=0.25
) +
ggtitle("PDAC Cell Type Atlas")


ggsave(
  "results/Figures/Figure1_UMAP_CellTypes.png",
  p1,
  width=12,
  height=9,
  dpi=600,
  bg="white"
)



############################################################
# Figure 2
# UMAP Samples
############################################################

p2 <- DimPlot(
  pdac,
  reduction="umap",
  group.by="sample",
  pt.size=0.25
)+
ggtitle("UMAP by Sample")


ggsave(
 "results/Figures/Figure2_UMAP_Samples.png",
 p2,
 width=10,
 height=8,
 dpi=600,
 bg="white"
)



############################################################
# Figure 3
# Tumor vs Normal
############################################################

p3 <- DimPlot(
 pdac,
 reduction="umap",
 group.by="condition",
 pt.size=0.25
)+
ggtitle("Tumor versus Normal")


ggsave(
 "results/Figures/Figure3_Tumor_Normal.png",
 p3,
 width=10,
 height=8,
 dpi=600,
 bg="white"
)



############################################################
# Marker Selection
############################################################

markers <- read.csv(
 "results/Tables/Cluster_Markers.csv"
)


clean_markers <- markers %>%

filter(
 p_val_adj < 0.01
) %>%

filter(
 !grepl(
 "^MT-|^RPS|^RPL|^LINC|^AC0|^RP11|^RNU|^LTB$",
 gene,
 ignore.case=TRUE
 )
) %>%

filter(
 pct.1 > 0.35,
 pct.2 < 0.20
) %>%

group_by(cluster) %>%

arrange(
 desc(avg_log2FC)
) %>%

slice_head(
 n=1
) %>%

ungroup()



write.csv(
 clean_markers,
 "results/Tables/Clean_Top_Markers_For_Plots.csv",
 row.names=FALSE
)



genes <- clean_markers$gene


genes <- genes[
 genes %in% rownames(pdac)
]


############################################################
# Scale only selected genes
############################################################

pdac <- ScaleData(
 pdac,
 features=genes,
 verbose=FALSE
)



############################################################
# Figure 4
# Marker Heatmap
############################################################


heat_pdac <- subset(
 pdac,
 downsample=50
)


heat <- DoHeatmap(
 heat_pdac,
 features=genes,
 group.by="celltype",
 size=4,
 angle=45
)+
NoLegend()+
scale_fill_gradient2(
 low="#2166ac",
 mid="white",
 high="#b2182b",
 midpoint=0
)+
ggtitle(
 "Top Cell Type Specific Marker Genes"
)+
theme(
 plot.title =
 element_text(
 hjust=0.5,
 size=16,
 face="bold",
 margin=margin(b=20)
 ),
 axis.text.y =
 element_text(
 size=10,
 face="italic"
 ),
 axis.text.x =
 element_blank(),
 axis.ticks.x =
 element_blank()
)



ggsave(
 "results/Figures/Figure4_Marker_Heatmap.png",
 heat,
 width=12,
 height=8,
 dpi=600,
 bg="white"
)



############################################################
# Figure 5
# Dot Plot
############################################################


dot <- DotPlot(
 pdac,
 features=genes,
 group.by="celltype",
 dot.scale=6
)+
RotatedAxis()+
ggtitle(
 "Marker Gene Expression Across Cell Types"
)+
theme(
 plot.title=
 element_text(
 hjust=0.5,
 size=16,
 face="bold"
 ),
 axis.text.x=
 element_text(
 angle=45,
 hjust=1,
 size=10,
 face="italic"
 )
)



ggsave(
 "results/Figures/Figure5_Marker_DotPlot.png",
 dot,
 width=12,
 height=8,
 dpi=600,
 bg="white"
)



############################################################
# Figure 6
# Cell Composition
############################################################


cellcount <- pdac@meta.data %>%
count(celltype)



p6 <- ggplot(
 cellcount,
 aes(
 x=reorder(celltype,n),
 y=n
 )
)+
geom_col(
 fill="steelblue"
)+
coord_flip()+
labs(
 title="Cell Type Composition",
 x="Cell Type",
 y="Number of Cells"
)



ggsave(
 "results/Figures/Figure6_Cell_Composition.png",
 p6,
 width=10,
 height=8,
 dpi=600,
 bg="white"
)



############################################################
# Figure 7
# QC Plot
############################################################


cells_keep <- sample(
 colnames(pdac),
 min(5000,ncol(pdac))
)


pdac_small <- subset(
 pdac,
 cells=cells_keep
)


qc <- VlnPlot(
 pdac_small,
 features=c(
 "nFeature_RNA",
 "nCount_RNA",
 "percent.mt"
 ),
 pt.size=0,
 ncol=3
)



ggsave(
 "results/Figures/Figure7_QC.png",
 qc,
 width=14,
 height=6,
 dpi=600,
 bg="white"
)



############################################################
# Figure 8
# Feature Plot
############################################################


feature_genes <- c(
 "EPCAM",
 "KRT19",
 "CD3D",
 "MS4A1",
 "LYZ",
 "COL1A1"
)


feature_genes <- feature_genes[
 feature_genes %in% rownames(pdac)
]


fp <- FeaturePlot(
 pdac,
 features=feature_genes,
 reduction="umap",
 order=TRUE,
 ncol=3
)



ggsave(
 "results/Figures/Figure8_FeaturePlot.png",
 fp,
 width=15,
 height=10,
 dpi=600,
 bg="white"
)



############################################################
# Save Visualization Object
############################################################


saveRDS(
 pdac,
 "data/processed/PDAC_Visualization.rds"
)


cat(
"
=====================================
Visualization Completed Successfully

Figures saved:
results/Figures/

=====================================
"
)
