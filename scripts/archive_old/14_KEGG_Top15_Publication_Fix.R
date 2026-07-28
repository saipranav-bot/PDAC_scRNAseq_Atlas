library(ggplot2)
library(dplyr)
library(stringr)


input <- "results/Pathways/KEGG_Tumor_Upregulated.csv"

output <- "results/Pathways/KEGG_Tumor_Top15_Publication_Fixed.png"


# Read KEGG enrichment
df <- read.csv(input)


# Select top 15 pathways
df <- df %>%
  arrange(p.adjust) %>%
  head(15)


# Convert GeneRatio
df$GeneRatio <- sapply(
  strsplit(df$GeneRatio, "/"),
  function(x){
    as.numeric(x[1]) / as.numeric(x[2])
  }
)


# Significance score
df$logP <- -log10(df$p.adjust)


# Wrap long pathway names
df$Description <- str_wrap(
  df$Description,
  width = 30
)


# Order
df$Description <- factor(
  df$Description,
  levels = rev(df$Description)
)



p <- ggplot(
  df,
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
  title = "Top 15 KEGG Pathways in PDAC Tumor Cells",
  x = "Gene Ratio",
  y = NULL
) +

theme_classic(
  base_size = 14
) +

theme(

  plot.title = element_text(
    hjust = 0.5,
    face = "bold",
    size = 17
  ),

  axis.text.y = element_text(
    size = 11,
    colour = "black"
  ),

  axis.text.x = element_text(
    size = 12
  ),

  legend.position = "right",

  plot.margin = margin(
    20,
    90,
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
  height = 8,
  dpi = 600,
  bg="white"
)


print(p)

cat("Saved:",output,"\n")
