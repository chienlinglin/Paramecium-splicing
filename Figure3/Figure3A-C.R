# Load library
library(tidyverse)
library(glue)

## 001. Figure 3A:  Alternative splicing events barplot
# Import rMATs summary data 
rMATs_summary <- read.csv("rmats_summary.csv")
# Plot
Fig3A <- gplot(transform(rMATs_summary,
                 group = factor(group, levels = c("AM11_g_w", "PM2_g_w", "PM5_g_w", "PM8_g_w", "PM11_g_w", "AM2_g_w", "AM5_g_w", "AM8_g_w"))), # Arrange the plot based on the timepoint order 
       aes(x = EventType, y = Value, fill = Compare, label = Value)) +
  geom_col() +
  geom_text(aes(label = ifelse(EventType %in% c("A3SS", "RI"), as.character(Value), "")), 
            position = position_stack(vjust = 0.5)) +
  facet_wrap(~ group, nrow = 2) +
  labs(
    x = "",
    y = "Number of AS events",
    fill = ""
  )+
  scale_fill_manual(values = c("#2A788EFF", "#FDE725FF"),
                    labels = c("green cells has higher PSI",
                               "white cells has higher PSI"
                    )) +
  theme_bw()

## 002. Figure 3C: DSIs barplot 
# Import DSI counts data 
sigDSI_counts <- read_csv("sigDSI_counts.csv")

# order the timepoint in plot
comparison_order <- c("AM11_g_w", "PM2_g_w", "PM5_g_w", "PM8_g_w", "PM11_g_w", "AM2_g_w", "AM5_g_w", "AM8_g_w") 
# Convert 'comparison' to a factor with the specified order
sigDSI_counts <- sigDSI_counts %>%
  mutate(comparison = factor(comparison, levels = comparison_order))

# Plot
Fig3C <- ggplot(sigDSI_counts, aes(x = comparison, y = count, fill = Difference)) +
  geom_col() +
  labs(x = "", y = "number of RIs", fill = "") +
  scale_fill_manual(values = c("#2A788EFF", "#FDE725FF", "#7AD151FF")) +
  geom_text(aes(label = as.character(count)), 
            position = position_stack(vjust = 0.5)) +
  theme(text = element_text(size = 14))
