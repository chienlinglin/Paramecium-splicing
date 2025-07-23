# Load library
library(tidyverse)
library(ggsignif)
library(ggpubr)

## 01. Figure 1A
# Import intron data 
Intron_data <- read_tsv("Intron_data.tsv")

# The middle point of intron = (end+start)/2
Intron_data <- Intron_data |>
  mutate(intron_midpoint=(intron_start+intron_end)/2,
         transcript_length = abs(three_prime_end-five_prime_end),
         intron_to_5_end = abs(intron_midpoint-five_prime_end)/transcript_length)

# Figure 1A
Fig1A <- ggplot(intron_position_data, aes(x=intron_to_5_end)) +
  geom_histogram(bins = 50,color="black", fill="lightblue")+
  labs(
    x = "Relative position to 5' end of gene",
    y = "Number of introns"
  )+
  theme_bw(base_size = 10)

## 02. Figure 1B: number of introns/ gene and gene expression level
# Import gene expression and intron number data
GE_intron_number <- read.csv("GE_intron_number.csv")
# Pivot longer
melt_GE_intron_number <- GE_intron_number |>
  pivot_longer(cols = 3:18, names_to = 'Sample', values_to = 'Normalized_counts')
# Plot
Fig1B <- ggplot(melt_GE_intron_number, aes(x=total_intron_number, y = log(Normalized_counts+1), fill = total_intron_number)) +
  geom_violin() +
  geom_boxplot(outliers = FALSE, width = 0.1)+ # outliners = FALSE: outliers are discarded and not included in the calculation
  geom_signif(
    comparisons = list(
      c("0", "1"),
      c("1", "2"),
      c("2", "3"),
      c("3", "4"),
      c("4", "5"),
      c("5", "6"),
      c("6", "7+")
    ),
    map_signif_level = TRUE,
    textsize = 6,
    margin_top = 0.08,
    step_increase = 0.02,
    tip_length = 0.01
  ) +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12))+
  theme(base_size = 35)+
  theme_bw()+
  labs(x = 'Number of introns per gene',
       y = 'Log2(normalized counts+1)',
       fill = "intron_number")

## 03. Figure 1C: First intron position and gene expression level 
# Import gene expression data 
avg_normalized_count <- read.csv("avg_normalized_count.csv")
# Find the nearest intron to 5 end in 1 gene
intron_position_min <- Intron_data |>
  group_by(Gene_ID) |>
  summarise(min_to_5_end = min(intron_to_5_end))

# Merge intron position with gene expression data 
GE_merge_5_end_distance <- intron_position_min |>
  left_join(avg_normalized_count, by = c("Gene_ID" = "gene_id"))

# Mutate the intron position
GE_merge_5_end_distance <- GE_merge_5_end_distance |>
  mutate(
    intron_position_category = case_when(
      min_to_5_end < 0.25 ~ "first_5_intron",
      min_to_5_end > 0.75 ~ "first_3_intron",
      TRUE ~ "middle"
    )
  )
# Melt the df 
melt_GE_5_end_distance <- GE_merge_5_end_distance|>
  pivot_longer(cols = 3:18, names_to = 'Sample', values_to = 'Normalized_counts')

# Plot
Fig1C <- ggplot(melt_GE_5_end_distance, aes(x=intron_position_category, y = log(Normalized_counts+1), fill = intron_position_category)) +
  geom_violin() +
  geom_boxplot(outliers = FALSE, width = 0.1)+ 
  geom_signif(
    comparisons = list(
      c("first_5_intron", "middle"),
      c("middle", "first_3_intron"), 
      c("first_5_intron", "first_3_intron")),
    map_signif_level = TRUE,
    textsize = 6,
    margin_top = 0.08,
    step_increase = 0.05,
    tip_length = 0.01
  ) +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12))+
  theme(base_size = 35)+ 
  theme_bw()+
  labs(x = 'Gene category based on first intron position',
       y = 'Log(normalized counts+1)',
       fill = "intron_position_category")
