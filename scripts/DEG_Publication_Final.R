############################################################
# DEG_Publication_Final.R
# Publication Quality DEG Analysis
# PDAC scRNA-seq Tumor vs Normal
############################################################


library(Seurat)
library(ggplot2)
library(dplyr)
library(ggrepel)
library(pheatmap)


set.seed(123)


dir.create(
"results/DEG",
recursive=TRUE,
showWarnings=FALSE
)


############################################################
# Load data
############################################################


pdac <- readRDS(
"data/processed/PDAC_Annotated.rds"
)


DefaultAssay(pdac)="RNA"


pdac[["RNA"]] <- JoinLayers(
pdac[["RNA"]]
)



############################################################
# DEG table
############################################################


deg <- read.csv(
"results/DEG/Tumor_vs_Normal_DEG.csv"
)


############################################################
# Volcano preparation
############################################################


minp <- min(
deg$p_val_adj[deg$p_val_adj>0],
na.rm=TRUE
)


deg <- deg %>%
mutate(

p_adj =
ifelse(
p_val_adj==0,
minp/10,
p_val_adj
),


Significance =
case_when(

p_val_adj <0.05 &
avg_log2FC >0.5
~ "Tumor",

p_val_adj <0.05 &
avg_log2FC < -0.5
~ "Normal",

TRUE
~ "NS"

),

neglog =
-log10(p_adj)

)



top_labels <- deg %>%
filter(Significance!="NS") %>%
arrange(p_val_adj) %>%
head(20)



############################################################
# FIGURE 1 Volcano
############################################################


volcano <- ggplot(
deg,
aes(
avg_log2FC,
neglog,
color=Significance
)
)+


geom_point(
size=2,
alpha=0.7
)+


scale_color_manual(
values=c(
"Tumor"="#d73027",
"Normal"="#4575b4",
"NS"="grey75"
)
)+


geom_vline(
xintercept=c(-0.5,0.5),
linetype="dashed"
)+


geom_hline(
yintercept=-log10(0.05),
linetype="dashed"
)+


geom_text_repel(
data=top_labels,
aes(label=gene),
size=4,
max.overlaps=20
)+


theme_classic(
base_size=16
)+


theme(
plot.title=
element_text(
hjust=.5,
face="bold"
)
)+


labs(
title="Differential Expression: Tumor vs Normal PDAC",
x="log2 Fold Change",
y="-log10 Adjusted P-value"
)



ggsave(
"results/DEG/Figure9_DEG_Volcano.png",
volcano,
width=10,
height=8,
dpi=600,
bg="white"
)




############################################################
# FIGURE 2 Heatmap
############################################################


tumor <- deg %>%
filter(Significance=="Tumor") %>%
arrange(desc(avg_log2FC)) %>%
head(10)


normal <- deg %>%
filter(Significance=="Normal") %>%
arrange(avg_log2FC) %>%
head(10)



top20 <- bind_rows(
tumor,
normal
)


write.csv(
top20,
"results/DEG/Top20_DEG_Selected.csv",
row.names=FALSE
)



genes <- top20$gene



expr <- AverageExpression(
pdac,
features=genes,
group.by="sample"
)$RNA



expr <- t(scale(t(expr)))



annotation <- data.frame(
Condition=
pdac@meta.data %>%
select(sample,condition) %>%
distinct() %>%
arrange(match(sample,colnames(expr))) %>%
pull(condition)
)


rownames(annotation)=colnames(expr)



png(
"results/DEG/Figure10_Top20_DEG_Heatmap.png",
3000,
2400,
res=300
)


pheatmap(

expr,

annotation_col=annotation,

cluster_cols=FALSE,

cluster_rows=TRUE,

color=
colorRampPalette(
c(
"#2166ac",
"white",
"#b2182b"
)
)(100),

fontsize_row=12,

fontsize_col=12,

border_color=NA,

main=
"Top Differentially Expressed Genes\nTumor vs Normal PDAC"

)


dev.off()



############################################################
# FIGURE 3 DEG DotPlot
############################################################


Idents(pdac)="condition"


dot <- DotPlot(
pdac,
features=genes,
group.by="condition"
)+

RotatedAxis()+

ggtitle(
"Expression of Top Differential Genes"
)+

theme_classic(
base_size=15
)


ggsave(

"results/DEG/Figure11_DEG_DotPlot.png",

dot,

width=10,

height=7,

dpi=600,

bg="white"

)



cat(
"
===================================
FINAL DEG FIGURES GENERATED
===================================

Figure9_DEG_Volcano.png

Figure10_Top20_DEG_Heatmap.png

Figure11_DEG_DotPlot.png

===================================
"
)
