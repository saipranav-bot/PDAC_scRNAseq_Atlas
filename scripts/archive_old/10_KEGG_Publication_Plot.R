library(ggplot2)
library(dplyr)
library(stringr)

# Input
input <- "results/Pathways/KEGG_Tumor_Upregulated.csv"

# Output
output <- "results/Pathways/KEGG_Tumor_Publication_Fixed.png"


# Load data
kegg <- read.csv(input)


# Select top 12 pathways (publication readable)
kegg_top <- kegg %>%
  arrange(p.adjust) %>%
  head(12)


# Calculate GeneRatio numeric
kegg_top$GeneRatio <- sapply(
  strsplit(kegg_top$GeneRatio, "/"),
  function(x){
    as.numeric(x[1])/as.numeric(x[2])
  }
)


# Convert p value
kegg_top$logP <- -log10(kegg_top$p.adjust)


# Wrap pathway names
kegg_top$Description <- str_wrap(
  kegg_top$Description,
  width = 32
)


# Order pathways
kegg_top$Description <- factor(
  kegg_top$Description,
  levels = rev(kegg_top$Description)
)



p <- ggplot(
  kegg_top,
  aes(
    x = GeneRatio,
    y = Description,
    size = Count,
    color = logP
  )
) +

geom_point(
  alpha = 0.85
) +


scale_color_gradient(
  low = "#2166AC",
  high = "#B2182B",
  name = "-log10\nadjusted p"
) +


scale_size_continuous(
  range = c(4,10),
  name = "Gene count"
) +


labs(
  title = "KEGG Pathway Enrichment in PDAC Tumor Cells",
  x = "Gene Ratio",
  y = NULL
) +


theme_classic(
  base_size = 15
) +


theme(

  plot.title = element_text(
    hjust = 0.5,
    face = "bold",
    size = 18
  ),

  axis.text.y = element_text(
    size = 12,
    colour = "black"
  ),

  axis.text.x = element_text(
    size = 12,
    colour = "black"
  ),

  legend.title = element_text(
    face="bold"
  ),

  legend.position="right",

  plot.margin = margin(
    20,
    80,
    20,
    20
  )
) +


coord_cartesian(
  clip="off"
)



ggsave(
  output,
  p,
  width = 11,
  height = 7,
  dpi = 600,
  bg="white"
)


print(p)

cat("Saved:",output,"\n")
