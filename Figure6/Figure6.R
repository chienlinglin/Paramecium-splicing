# Load library
library(tidyverse)
library(ggsignif)

## 001. Figure 6B. Intron age and PSI
# Import intron PSI and intron age data
Pb_intron_age <- read_csv("PSI_intron_age.csv")

# Import intron PSI data
PSI_data <- read_csv("PSI_data.csv")

# Import intron id to gene id map
intron_id_to_gene_id_all <- read_csv("intron_id_to_gene_id_all.csv")

# Merge intron PSI data and intron age data
Pb_intron_age_merge_PSI <- PSI_data |>
  inner_join(Pb_intron_age, by = c("intron_id", "GeneID" = "gene_id"))

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
# Plot 
Pb_intron_age$intron_age <- as.factor(Pb_intron_age$intron_age)
Fig6C <- ggplot(Pb_intron_age, aes(x = as.factor(intron_age), y = intron_gc_content, fill = as.factor(intron_age))) + 
  geom_boxplot(outliers = FALSE, width = 0.3) +
  labs(
    y = "intron content (%)",
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
      c("1", "3")
    ),
    map_signif_level = TRUE,
    textsize = 6,
    margin_top = -0.35,
    step_increase = 0.05,
    tip_length = 0.01
  )

## 003. Figure 6D. Intron age and intron length
# Perform KS test bw 1 & 2, 2 & 3 intron age group length
# Convert intron_age to factor
Pb_intron_age$intron_age <- as.factor(Pb_intron_age$intron_age)

ks_results <- list(
  "1_vs_2" = ks.test(Pb_intron_age$intron_width[Pb_intron_age$intron_age == "1"],
                     Pb_intron_age$intron_width[Pb_intron_age$intron_age == "2"], alternative = "two.sided", exact = FALSE),
  "1_vs_3" = ks.test(Pb_intron_age$intron_width[Pb_intron_age$intron_age == "1"],
                     Pb_intron_age$intron_width[Pb_intron_age$intron_age == "3"], alternative = "two.sided", exact = FALSE)
)

# Extract p-values
ks_p_values <- data.frame(
  comparison = c("1_vs_2", "1_vs_3"),
  p_value = c(ks_results$`1_vs_2`$p.value, ks_results$`1_vs_3`$p.value)
)

# Plot
Pb_intron_age$intron_age <- as.factor(Pb_intron_age$intron_age)
Fig6D <- ggplot(Pb_intron_age, aes(x = as.factor(intron_age), y = intron_width, fill = as.factor(intron_age))) + 
  geom_boxplot(outliers = FALSE, width = 0.2) +
  labs(
    y = "intron length (bp)",
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
      c("1", "3")
    ),
    annotations = sprintf("p = %.3g", ks_p_values$p_value),  # Use KS test p-values
    textsize = 6,
    margin_top = -0.75,
    step_increase = 0.02,
    tip_length = 0.01
  ) 

## 004. Figure 6E. Maximum intron age and gene expression level
# Import gene expression data 
avg_normalized_count <- read_csv('avg_normalized_count.csv')

# Summarize the oldest intron age in each gene 
Pb_intron_age_gene <- Pb_intron_age |>
  group_by(gene_id) |>
  summarize(max_intron_age = max(as.numeric(intron_age), na.rm = TRUE))

# merge intron age data and gene expression 
Pb_intron_age_merge_GE <- Pb_intron_age_gene |>
  inner_join(avg_normalized_count, by = "gene_id")

Pb_intron_age_merge_GE_long <- Pb_intron_age_merge_GE |>
  pivot_longer(
    cols = starts_with("avg"), 
    names_to = "samples", 
    values_to = "avg_GE"
  )

Pb_intron_age_merge_GE_long$max_intron_age <- as.factor(Pb_intron_age_merge_GE_long$max_intron_age)

# Plot 
Fig6E <- ggplot(Pb_intron_age_merge_GE_long, aes(x = as.factor(max_intron_age), y = log(avg_GE+1), fill = as.factor(max_intron_age))) + 
  geom_violin()+
  geom_boxplot(outliers = FALSE, width = 0.1) +
  labs(
    y = "log(Normalized count +1)",
    x = "max_intron_age", 
    fill = 'max_intron_age'
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
    margin_top = 0.08,
    step_increase = 0.02,
    tip_length = 0.01
  )