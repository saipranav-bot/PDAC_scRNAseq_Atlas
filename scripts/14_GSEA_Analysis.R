library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)


dir.create(
"results/Pathways/GSEA",
recursive = TRUE,
showWarnings = FALSE
)


# Load DEG table

deg <- read.csv(
"results/DEG/Significant_Tumor_vs_Normal_DEG.csv"
)


# Create ranked gene list

gene_rank <- deg %>%
  select(gene, avg_log2FC) %>%
  filter(!is.na(avg_log2FC))


# Convert symbols

gene_id <- bitr(
gene_rank$gene,
fromType="SYMBOL",
toType="ENTREZID",
OrgDb=org.Hs.eg.db
)


ranked <- gene_rank %>%
  inner_join(
    gene_id,
    by=c("gene"="SYMBOL")
  )


geneList <- ranked$avg_log2FC

names(geneList) <- ranked$ENTREZID


geneList <- sort(
geneList,
decreasing = TRUE
)



################################################
# GO GSEA
################################################


gsea_go <- gseGO(
geneList = geneList,
OrgDb = org.Hs.eg.db,
ont="BP",
keyType="ENTREZID",
minGSSize=10,
maxGSSize=500,
pvalueCutoff=0.05,
verbose=FALSE
)


write.csv(
as.data.frame(gsea_go),
"results/Pathways/GSEA/GSEA_GO_Tumor.csv",
row.names=FALSE
)



png(
"results/Pathways/GSEA/GSEA_GO_Dotplot.png",
width=3000,
height=2500,
res=300
)

print(
dotplot(
gsea_go,
showCategory=15
)+
ggtitle(
"GO GSEA - PDAC Tumor vs Normal"
)
)

dev.off()



################################################
# KEGG GSEA
################################################


gsea_kegg <- gseKEGG(
geneList=geneList,
organism="hsa",
minGSSize=10,
maxGSSize=500,
pvalueCutoff=0.05,
verbose=FALSE
)


write.csv(
as.data.frame(gsea_kegg),
"results/Pathways/GSEA/GSEA_KEGG_Tumor.csv",
row.names=FALSE
)



png(
"results/Pathways/GSEA/GSEA_KEGG_Dotplot.png",
width=3000,
height=2500,
res=300
)


print(
dotplot(
gsea_kegg,
showCategory=15
)+
ggtitle(
"KEGG GSEA - PDAC Tumor vs Normal"
)
)


dev.off()



cat("GSEA completed\n")
