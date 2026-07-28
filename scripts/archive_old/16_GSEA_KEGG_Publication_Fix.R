library(ggplot2)
library(dplyr)
library(stringr)


input <- "results/Pathways/GSEA/GSEA_KEGG_Tumor.csv"

output <- "results/Pathways/GSEA/GSEA_KEGG_Publication_Fixed.png"


# Load data
df <- read.csv(input)


# Select top pathways
df <- df %>%
  arrange(p.adjust) %>%
  head(10)


# Calculate significance
df$logP <- -log10(df$p.adjust)


# Wrap names
df$Description <- str_wrap(
  df$Description,
  width = 30
)


# Order pathways
df$Description <- factor(
  df$Description,
  levels = rev(df$Description)
)



p <- ggplot(
  df,
  aes(
    x = NES,
    y = Description,
    size = setSize,
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
  range = c(4,9),
  name = "Gene set size"
) +


labs(
  title = "KEGG GSEA Enrichment - PDAC Tumor vs Normal",
  x = "Normalized Enrichment Score (NES)",
  y = NULL
) +


theme_classic(
  base_size = 14
) +


theme(

  plot.title = element_text(
    hjust = 0.5,
    face = "bold",
    size = 16
  ),

  axis.text.y = element_text(
    size = 10,
    colour = "black"
  ),

  axis.text.x = element_text(
    size = 11,
    colour = "black"
  ),

  legend.title = element_text(
    face = "bold"
  ),

  legend.text = element_text(
    size = 10
  ),

  plot.margin = margin(
    20,
    100,
    20,
    40
  )
) +


coord_cartesian(
  clip = "off"
)



ggsave(
  output,
  p,
  width = 12,
  height = 7,
  dpi = 600,
  bg = "white"
)


print(p)

cat("Saved:", output, "\n")
