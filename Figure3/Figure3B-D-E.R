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
# Import A3SS data
merged_A3SS_data <- read_csv("erged_A3SS_data.csv")

# Calculate the distance bw canonical 3'SS and A3SS 3'SS
merged_A3SS_data <- merged_A3SS_data |>
  mutate(distance_canSS_ASS = ifelse(strand == "+", abs(can_exon_start - A3_exon_start), abs(can_exon_end - A3_exon_end)))

A3SS_distance <- merged_A3SS_data$distance_canSS_ASS
# Create the cut variable with all levels explicitly
breaks <- c(seq(0, 15, 1), Inf)
labels <- c(as.character(0:14), '>14')
A3SS_distance.cut <- cut(A3SS_distance, breaks = breaks, labels = labels, right = FALSE)

# Force all levels to be present
A3SS_distance.cut <- factor(A3SS_distance.cut, levels = labels)

# Convert to a df for ggplot
data <- data.frame(A3SS_distance.cut)

# Plot using ggplot2
Fig3D <- ggplot(data, aes(x = A3SS_distance.cut)) +
  geom_bar(fill = 'black', color = 'black') +
  scale_x_discrete(breaks = labels,
                   labels = labels) +
  labs(x = 'Intron width', y = 'Counts', title = '') +
  theme(axis.text.x = element_text(color = "black", size = 9),
        axis.text.y = element_text(color = "black", size = 9))+
  theme(base_size = 32)+
  scale_x_discrete(drop = FALSE) +
  theme_bw()