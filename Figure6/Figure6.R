# Load library
library(tidyverse)
library(ggsignif)

## 001. Figure 6B. Intron age and PSI
# Import intron PSI and intron age data
Intron_age <- read_csv("PSI_intron_age.csv")

# Import intron PSI data
PSI_data <- read_csv("PSI_data.csv")

# Import intron id to gene id map
intron_id_to_gene_id_all <- read_csv("intron_id_to_gene_id_all.csv")

# Merge intron PSI data and intron age data
Pb_intron_age_merge_PSI <- PSI_data |>
  inner_join(Intron_age, by = c("intron_id", "GeneID" = "gene_id"))

# melt intron age merged PSI df
Pb_intron_age_merge_PSI_long <- Pb_intron_age_merge_PSI |>
  pivot_longer(
    cols = starts_with("avg"), 
    names_to = "samples", 
    values_to = "avg_PSI"
  )

# plot 
Fig6B <- ggplot(Pb_intron_age_merge_PSI, aes(x = as.factor(intron_age), y = avg_PSI*100, fill = as.factor(intron_age))) + 
  geom_boxplot(outliers = FALSE, width = 0.2) +
  labs(
    y = "PSI (%)",
    x = "intron_age", 
    fill = 'intron_age'
  ) +
  theme_bw(base_size = 12) +
  guides(fill = guide_legend(title = "")) +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12)) +
  geom_signif(
    comparisons = list(
      c("1", "2"),
      c("2", "3")
    ),
    map_signif_level = TRUE,
    textsize = 6,
    margin_top = -0.85,
    step_increase = 0.02,
    tip_length = 0.01
  ) 
## 002. Figure 6C. Intron age and Intron GC content

## 003. Figure 6D. Intron age and intron length

## 004. Figure 6E. Maximum intron age and gene expression level
