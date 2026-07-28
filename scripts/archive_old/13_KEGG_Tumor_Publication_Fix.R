library(ggplot2)
library(dplyr)
library(stringr)

input <- "results/Pathways/KEGG_Tumor_Upregulated.csv"

output <- "results/Pathways/KEGG_Tumor_Publication_Fixed.png"


# Load data
kegg <- read.csv(input)


# Select top 10 pathways
kegg_top <- kegg %>%
  arrange(p.adjust) %>%
  head(10)


# Convert pathway order
kegg_top$Description <- factor(
  kegg_top$Description,
  levels = rev(kegg_top$Description)
)


# Wrap names
kegg_top$Description <- str_wrap(
  kegg_top$Description,
  width = 30
)


# -log10 adjusted p
kegg_top$logP <- -log10(kegg_top$p.adjust)



p <- ggplot(
  kegg_top,
  aes(
    x = logP,
    y = Description,
    size = Count
  )
) +

geom_point(
  color = "#2166AC",
  alpha = 0.85
) +

scale_size_continuous(
  range=c(4,10),
  name="Gene Count"
) +

labs(
  title="KEGG Pathway Enrichment in PDAC Tumor Cells",
  x=expression(-log[10](adjusted~p)),
  y=NULL
) +

theme_classic(base_size=14) +

theme(

  plot.title = element_text(
    hjust=0.5,
    face="bold",
    size=16
  ),

  axis.text.y = element_text(
    size=12,
    colour="black"
  ),

  axis.text.x = element_text(
    size=12,
    colour="black"
  ),

  legend.title = element_text(
    face="bold"
  ),

  plot.margin = margin(
    15,40,15,120
  )
)



ggsave(
 output,
 p,
 width=9,
 height=7,
 dpi=600
)


cat("Saved:",output,"\n")
