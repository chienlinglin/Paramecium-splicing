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
