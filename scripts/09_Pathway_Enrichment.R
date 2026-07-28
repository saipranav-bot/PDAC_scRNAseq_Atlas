############################################################
# 09_Pathway_Enrichment.R
# GO + KEGG Functional Enrichment of Tumor vs Normal DEGs
# PDAC scRNA-seq
############################################################

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(dplyr)
library(tibble)

set.seed(123)

############################################################
# Create output folder
############################################################

dir.create("results/Pathways", recursive = TRUE, showWarnings = FALSE)

############################################################
# Load DEG results
############################################################

deg <- read.csv(
  "results/DEG/Significant_Tumor_vs_Normal_DEG.csv"
)

############################################################
# Separate Tumor and Normal enriched genes
############################################################

tumor_genes <- deg %>%
  filter(avg_log2FC > 0.5,
         p_val_adj < 0.05) %>%
  pull(gene)


normal_genes <- deg %>%
  filter(avg_log2FC < -0.5,
         p_val_adj < 0.05) %>%
  pull(gene)


cat("Tumor enriched genes:", length(tumor_genes), "\n")
cat("Normal enriched genes:", length(normal_genes), "\n")


############################################################
# Convert gene symbols to ENTREZ IDs
############################################################

convert_genes <- function(gene_list){

  bitr(
    gene_list,
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  ) %>%
    distinct(SYMBOL, .keep_all = TRUE)

}


tumor_id <- convert_genes(tumor_genes)
normal_id <- convert_genes(normal_genes)


############################################################
# Save converted genes
############################################################

write.csv(
  tumor_id,
  "results/Pathways/Tumor_Genes_Entrez.csv",
  row.names = FALSE
)

write.csv(
  normal_id,
  "results/Pathways/Normal_Genes_Entrez.csv",
  row.names = FALSE
)


############################################################
# Function for GO + KEGG analysis
############################################################

run_enrichment <- function(entrez,
                           name){

###########################################################
# GO Biological Process
###########################################################

ego <- enrichGO(
  gene = entrez,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)


write.csv(
  as.data.frame(ego),
  paste0(
    "results/Pathways/GO_",
    name,
    ".csv"
  ),
  row.names = FALSE
)


###########################################################
# KEGG pathway
###########################################################

ekegg <- enrichKEGG(
  gene = entrez,
  organism = "hsa",
  pvalueCutoff = 0.05
)


ekegg <- setReadable(
  ekegg,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID"
)


write.csv(
  as.data.frame(ekegg),
  paste0(
    "results/Pathways/KEGG_",
    name,
    ".csv"
  ),
  row.names = FALSE
)



###########################################################
# GO Dotplot
###########################################################

if(nrow(as.data.frame(ego)) > 0){

png(
paste0(
"results/Pathways/GO_",
name,
"_Dotplot.png"
),
width = 3000,
height = 2200,
res = 300
)

print(
dotplot(
ego,
showCategory = 15,
font.size = 12
)+
ggtitle(
paste(
"GO Biological Processes -",
name
)
)
)

dev.off()

}


###########################################################
# KEGG Dotplot
###########################################################

if(nrow(as.data.frame(ekegg)) > 0){

png(
paste0(
"results/Pathways/KEGG_",
name,
"_Dotplot.png"
),
width = 3000,
height = 2200,
res = 300
)


print(
dotplot(
ekegg,
showCategory = 15,
font.size = 12
)+
ggtitle(
paste(
"KEGG Pathways -",
name
)
)
)


dev.off()

}


}


############################################################
# Run enrichment
############################################################


run_enrichment(
tumor_id$ENTREZID,
"Tumor_Upregulated"
)


run_enrichment(
normal_id$ENTREZID,
"Normal_Upregulated"
)



cat("
========================================
Pathway enrichment completed successfully

Outputs saved:
results/Pathways/

========================================
")
