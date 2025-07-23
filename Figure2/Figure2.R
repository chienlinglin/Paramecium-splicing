# Load library
library(tidyverse)
library(RColorBrewer)
library(VennDiagram)
library(ggVennDiagram)
library(ggvenn)

## 001. Figure 2A: Upset plot
# Import RBHs data 
spliceosome_data_matrix <- read.csv("spliceosome_data_matrix")

# Create a list of sets for each species
sets_list <- list(
  "H. sapiens" = spliceosome_data_matrix$`Gene Symbol`[spliceosome_data_matrix$`H. sapiens` == 1],
  "M. musculus" = spliceosome_data_matrix$`Gene Symbol`[spliceosome_data_matrix$`Mus_musculus.fa` == 1],
  "D. rerio" = spliceosome_data_matrix$`Gene Symbol`[spliceosome_data_matrix$`Danio_rerio.fa` == 1],
  "D. melanogaster" = spliceosome_data_matrix$`Gene Symbol`[spliceosome_data_matrix$`Drosophila_melanogaster.fa` == 1],
  "C. elegans" = spliceosome_data_matrix$`Gene Symbol`[spliceosome_data_matrix$`C_elegans.fa` == 1],
  "S. cerevisiae" = spliceosome_data_matrix$`Gene Symbol`[spliceosome_data_matrix$`Saccharomyces_cerevisiae.fa` == 1],
  "P. bursaria" = spliceosome_data_matrix$`Gene Symbol`[spliceosome_data_matrix$`Paramecium_bursaria_DK2_uniq_protein.faa` == 1],
  "P. tetraurelia" = spliceosome_data_matrix$`Gene Symbol`[spliceosome_data_matrix$`Paramecium_tetraurelia.fa` == 1],
  "T. thermophila" = spliceosome_data_matrix$`Gene Symbol`[spliceosome_data_matrix$`Tetrahymena_thermophila.fa` == 1],
  "A. thaliana" = spliceosome_data_matrix$`Gene Symbol`[spliceosome_data_matrix$`Arabidopsis_thaliana.fa` == 1]
)

# Plot
Fig2A <- ggVennDiagram(sets_list, order.set.by = "none")

## 002. Figure 2B

