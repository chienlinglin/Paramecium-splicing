# Load library
library(tidyverse)
library(plotrix)
library(ggpubr)

## 001. Figure 4A: GC content and PSI correlation
# Import PSI - GC content
PSI_gc_content <- read_csv("PSI_gc_content.csv")

# Round GC content to one decimal place
PSI_gc_content$GC_content_round <- round((PSI_gc_content$intron_gc_content)/100, 1)

# pivor longer 
PSI_gc_content_sum <- PSI_gc_content |> 
  pivot_longer(
    cols = starts_with("avg_PSI"), 
    names_to = "samples", 
    values_to = "avg_PSI"
  )

# Summary the PSI in each GC content group 
mean_PSI_sum <- PSI_gc_content_sum %>%
  group_by(GC_content_round) %>%
  dplyr::summarize(mean_PSI = mean(avg_PSI, na.rm = TRUE),
                   se_PSI = std.error(avg_PSI, na.rm = TRUE))

# Plot
Fig4A <- ggplot(mean_PSI_sum, aes(x = GC_content_round)) +
  geom_line(aes(y = mean_PSI, color = "mean PSI"), linewidth = 1.5) +
  geom_errorbar(aes(y = mean_PSI, ymin = mean_PSI - se_PSI, ymax = mean_PSI + se_PSI), color = 'black', width = 0.015) +
  labs(x = "GC content", y = "PSI", color = "Sample")

## 002. Intron retention and intron length 
# Import intron length and PSI data 
intron_length_PSI <- read_csv("intron_length_PSI.csv")

# melt data 
intron_length_PSI_all_long <- intron_length_PSI_all |>
  pivot_longer(
    cols = starts_with("avg"),
    names_to = "samples",
    values_to = "avg_PSI"
  )

# Divide intron length into: 15-22, 23-26, >=27
intron_length_PSI_all_long <- intron_length_PSI_all_long |>
  mutate(intron_length_cat = case_when(intron_width <23 ~ "15-22",
                                       intron_width < 27 ~ "23-26",
                                       intron_width >=27 ~ ">=27"))

# Order factor level based on the length 
intron_length_PSI_all_long$intron_length_cat <- factor(
  intron_length_PSI_all_long$intron_length_cat,
  levels = c("15-22", "23-26", ">=27"),
  labels = c("15-22", "23-26", ">=27")
)

# Change >= symbol
levels(intron_length_PSI_all_long$intron_length_cat)[3] <- "≥27"
# Plot
Fig4B <- ggplot(intron_length_PSI_all_long, aes(x = intron_length_cat, y = avg_PSI, fill = intron_length_cat)) +
  geom_boxplot(outliers = FALSE, width = 0.2) +
  #stat_n_text(angle = 45) + 
  geom_signif(
    comparisons = list(
      c("15-22", "23-26"),
      c("23-26", "≥27")
    ),
    map_signif_level = TRUE,
    textsize = 6,
    margin_top = -0.8,
    step_increase = 0.02,
    tip_length = 0.01
  ) +
  theme_bw(base_size = 12)

## 003. Figure 4C: correlation bw PSI and gene expression level
# Import PSI data: avg PSI in gene
avg_PSI_each_gene <- read_csv("avg_PSI_each_gene.csv")

# Import gene expression data
avg_normalized_counts <- read_csv("avg_normalized_counts.csv")

# Merge the gene expression data to PSI 
merge_GE_PSI <- avg_normalized_count |>
  right_join(average_PSI_per_gene, by = "gene_id")

# melt the df: pivot longer for GE columns and PSI columns 
merge_GE_PSI_long <- merge_GE_PSI |>
  pivot_longer(
    cols = starts_with("avg_"), 
    names_to = "sample_GE", 
    values_to = "normalized_counts"
  ) |>
  pivot_longer(
    cols = starts_with("mean_PSI"), 
    names_to = "sample_PSI", 
    values_to = "mean_PSI"
  )

merge_GE_PSI_long_1$log_PSI <- log(merge_GE_PSI_long_1$mean_PSI + 0.01)
merge_GE_PSI_long_1$log_GE <- log(merge_GE_PSI_long_1$normalized_counts + 1)

# Plot
Fig4C <- ggscatter(
  merge_GE_PSI_long_1,
  x = "log_PSI",
  y = "log_GE",
  add = "reg.line",
  conf.int = FALSE,
  palette = "#A1E3F9",           
  add.params = list(color = "#8B0000"),
  cor.coef = TRUE,
  cor.method = "pearson",
  alpha = 1/200,
  size = 0.0015
) +     
ggsave("scatter_plot.png", p, width = 4, height = 4, units = "in", dpi = 600)

## 004. Figure 4D: intron position and PSI association
# Import intron position data 
intron_position_data <- read_tsv("/data/Data/Giang/Paramecium/Pb_intron_features/Intron_data.tsv")

# Select intron position column 
intron_position_data <- intron_position_data |>
  select(Gene_ID, seqnames, intron_start, intron_end, intron_to_5_end)

# Mutate intron position data to divide introns based on its position in gene
intron_position_data <- intron_position_data |>
  mutate(
    intron_position_category = case_when(
      intron_to_5_end <= 0.25 ~ "5_intron",
      intron_to_5_end >= 0.75 ~ "3_intron",
      TRUE ~ "middle_intron"
    )
  )

# Import intron PSI data
PSI_data <- read_csv("PSI_data.csv")
# Merge with PSI data 
intron_position_data <- intron_position_data |>
  mutate(intron_id = paste0(seqnames, "_", intron_start, "_", intron_end))
intron_position_merge_PSI <- intron_position_data |>
  inner_join(filtered_merged_data_1, by = c("intron_id"))
intron_position_merge_PSI <- intron_position_merge_PSI |>
  dplyr::select(intron_id, intron_to_5_end, intron_position_category, starts_with("avg"))
# melt PSI columns in intron position merge PSI df
intron_position_merge_PSI_long <- intron_position_merge_PSI |>
  pivot_longer(
    cols = starts_with("avg"), 
    names_to = "samples", 
    values_to = "avg_PSI"
  )

# Relevel
intron_position_merge_PSI_long$intron_position_category <- factor(
  intron_position_merge_PSI_long$intron_position_category,
  levels = c("5_intron", "middle_intron", "3_intron")
)

# Plot 
Fig4D <- ggplot(intron_position_merge_PSI_long, aes(x = intron_position_category, y = avg_PSI * 100, fill = intron_position_category)) + 
  geom_boxplot(outliers = FALSE, width = 0.2) +
  labs(
    y = "PSI (%)",
    x = "intron_position_category", 
    fill = 'intron_position_category'
  )+
  theme_bw(base_size = 12) +
  guides(fill = guide_legend(title = "")) +
  theme(axis.text.x = element_text(color = "black", size = 12),
        axis.text.y = element_text(color = "black", size = 12)) +
  geom_signif(
    comparisons = list(
      c("5_intron", "middle_intron"),
      c("middle_intron", "3_intron"),
      c("5_intron", "3_intron")
    ),
    map_signif_level = TRUE,
    textsize = 6,
    margin_top = -0.85,
    step_increase = 0.02,
    tip_length = 0.01
  )