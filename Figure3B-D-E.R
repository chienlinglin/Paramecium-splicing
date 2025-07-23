# Load library
library(tidyverse)

## 001. Figure 3B: Summary Sig A3SS events
# Import A3SS counts data 
sigA3SS_counts <- read_csv("")

# Order of sample on barplot
comparison_order <- c("AM11", "PM2", "PM5", "PM8",
                   "PM11", "AM2", "AM5", "AM8")

# Convert 'comparison' column to a factor with the specified order
sigA3SS_counts$timepoints <- factor(sigA3SS_counts$timepoints, levels = comparison_order)
ggplot(group_counts, aes(x = timepoints, y = count, fill = Difference)) +
  geom_col() +
  labs(x = "Comparison", y = "number of Sig A3SS events", fill = "") +
  scale_fill_manual(values = c("#2A788EFF", "#FDE725FF", "#7AD151FF")) +
  geom_text(aes(label = as.character(count)),
            position = position_stack(vjust = 0.5)) +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

## 002. Fig3D: distance between A3SS 3'SS and canonical 3'SS 
