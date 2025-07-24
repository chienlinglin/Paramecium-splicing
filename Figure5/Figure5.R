# Load library
library(tidyverse)
library(ggfortify)

## 001. Figure 5A. PCA plot 
# Import PSI data of DSIs for all replicates 
PSI_data_all_replicates_scale_sig <- read_csv("PSI_data_all_replicates_scale_sig.csv")

# Run PCA
pca_res <- prcomp(PSI_data_all_replicates_scale_sig)

# Add cell type and timepoints data 
PSI_data_all_replicates_scale_sig_1 <- as.data.frame(PSI_data_all_replicates_scale_sig) |>
  rownames_to_column(var = "sample") |>
  mutate(cell_type = ifelse(grepl("g", sample), "green", "white"), 
         timepoints = sapply(strsplit(sample, "_"), `[`, 3)) 

# Plot
Fig5A <- autoplot(pca_res, data = PSI_data_all_replicates_scale_sig, 
              colour = "cell_type", shape = "timepoints", size = 4,
              palette = c("white" = "#393E46", "green" = "#537D5D")) +
  theme_bw() +
  scale_shape_manual(values = c(16, 17, 15, 18, 8, 3, 4, 0)) +
  scale_colour_manual(values = c("white" = "#393E46", "green" = "#537D5D"))