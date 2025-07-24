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


## 002. Figure 6C. Intron age and Intron GC content

## 003. Figure 6D. Intron age and intron length

## 004. Figure 6E. Maximum intron age and gene expression level
