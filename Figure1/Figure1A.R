# Load library
library(tidyverse)

## 01. Figure 1A
# Import intron data 
Intron_data <- read_tsv("Intron_data.tsv")

# The middle point of intron = (end+start)/2
Intron_data <- Intron_data |>
  mutate(intron_midpoint=(intron_start+intron_end)/2,
         transcript_length = abs(three_prime_end-five_prime_end),
         intron_to_5_end = abs(intron_midpoint-five_prime_end)/transcript_length)

# Figure 1A
ggplot(intron_position_data, aes(x=intron_to_5_end)) +
  geom_histogram(bins = 50,color="black", fill="lightblue")+
  labs(
    x = "Relative position to 5' end of gene",
    y = "Number of introns"
  )+
  theme_bw(base_size = 10)

## 02. Figure 1B: number of introns/ gene and 
