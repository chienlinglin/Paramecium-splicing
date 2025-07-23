# Load library
library(tidyverse)
library(plotrix)

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

## 02. Intron retention and intron length 
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