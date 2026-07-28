library(ggplot2)
library(dplyr)

input <- "results/Pathways/KEGG_Tumor_Upregulated.csv"

output <- "results/Pathways/KEGG_Tumor_Plot6_Fixed.png"


df <- read.csv(input)


# Select only top pathways
df <- df %>%
  arrange(p.adjust) %>%
  head(8)


# Convert GeneRatio
df$GeneRatio <- sapply(
  strsplit(df$GeneRatio, "/"),
  function(x){
    as.numeric(x[1]) / as.numeric(x[2])
  }
)


# Arrange order
df$Description <- factor(
  df$Description,
  levels = rev(df$Description)
)



p <- ggplot(
  df,
  aes(
    x = GeneRatio,
    y = Description
  )
) +

geom_col(
  fill = "steelblue",
  width = 0.7
) +

geom_text(
  aes(label = Count),
  hjust = -0.2,
  size = 4
) +

labs(
  title = "Top KEGG Enriched Pathways in PDAC Tumor",
  x = "Gene Ratio",
  y = NULL
) +

theme_classic(base_size = 14) +

theme(
  plot.title = element_text(
    hjust = 0.5,
    face = "bold",
    size = 16
  ),

  axis.text.y = element_text(
    size = 11
  ),

  plot.margin = margin(10,40,10,10)
)


ggsave(
  filename = output,
  plot = p,
  width = 9,
  height = 5.5,
  dpi = 600
)


cat("Saved:", output, "\n")
