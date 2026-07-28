library(clusterProfiler)
library(ggplot2)
library(dplyr)

input <- "results/Pathways/KEGG_Tumor_Upregulated.csv"

out <- "results/Pathways/KEGG_Tumor_Plot5_Fixed.png"

df <- read.csv(input)

# select top pathways
df <- df %>%
  arrange(p.adjust) %>%
  head(10)

# convert ratio
df$GeneRatio <- sapply(strsplit(df$GeneRatio, "/"),
                       function(x) as.numeric(x[1])/as.numeric(x[2]))

df$Description <- factor(df$Description,
                         levels = rev(df$Description))


p <- ggplot(df,
            aes(x = GeneRatio,
                y = Description,
                size = Count,
                color = p.adjust)) +

  geom_point() +

  scale_color_gradient(
    low="red",
    high="blue",
    name="Adjusted p-value"
  ) +

  scale_size_continuous(
    name="Gene count"
  ) +

  labs(
    title="Top KEGG Pathways in PDAC Tumor Cells",
    x="Gene Ratio",
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
      size=11
    ),

    legend.position="right"
  )


ggsave(
  out,
  p,
  width=9,
  height=6,
  dpi=600
)

cat("Saved:",out,"\n")
