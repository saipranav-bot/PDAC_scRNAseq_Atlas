library(ggplot2)
library(dplyr)
library(stringr)


input <- "results/Pathways/GSEA/GSEA_GO_Tumor.csv"

output <- "results/Pathways/GSEA/GSEA_GO_Publication_Fixed.png"


df <- read.csv(input)


# Select significant pathways
df <- df %>%
  arrange(p.adjust) %>%
  head(10)


# Calculate plotting values
df$logP <- -log10(df$p.adjust)


# Wrap pathway names
df$Description <- str_wrap(
  df$Description,
  width = 35
)


# Order
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
  range=c(5,12),
  name="Gene set size"
) +

labs(
  title="GO GSEA Enrichment - PDAC Tumor vs Normal",
  x="Normalized Enrichment Score (NES)",
  y=NULL
) +

theme_classic(
  base_size=15
) +

theme(

plot.title = element_text(
  hjust=0.5,
  face="bold",
  size=18
),

axis.text.y = element_text(
 size=12,
 colour="black"
),

axis.text.x = element_text(
 size=12
),

legend.title = element_text(
 face="bold"
),

plot.margin = margin(
20,60,20,20
)

)



ggsave(
output,
p,
width=10,
height=6,
dpi=600,
bg="white"
)


print(p)

cat("Saved:",output,"\n")
