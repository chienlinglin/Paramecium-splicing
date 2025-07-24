# Load library
library(tidyverse)
library(pheatmap)

## 001. Figure 5E: beta coefficient of SF DEGs
# Import z scaled PSI of 992 DSIs
PSI_data_all_replicates_scale_sig <- read_csv("PSI_sigRIs_all_samples_z_scaled.csv")

# Import 
