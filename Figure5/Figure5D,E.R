# Load library
library(tidyverse)
library(pheatmap)
library(ggpubr)

## 001. Figure 5E: beta coefficient of SF DEGs
# Import z scaled PSI of 992 DSIs
PSI_data_all_replicates_scale_sig <- read_csv("PSI_sigRIs_all_samples_z_scaled.csv")

# Import z scaled of log SF DEGs
Spliceosomal_DEGs_log_z_scaled <- read_csv("Spliceosomal_DEGs_log_z_scaled.csv")

# empty matrix to store coefficient and p values of lm 
beta_matrix <- matrix(NA, nrow = ncol(PSI_data_all_replicates_scale_sig), ncol = ncol(Spliceosomal_DEGs_log_z_scaled))
beta1_matrix <- matrix(NA, nrow = ncol(PSI_data_all_replicates_scale_sig), ncol = ncol(Spliceosomal_DEGs_log_z_scaled))
p_beta_matrix <- matrix(NA, nrow = ncol(PSI_data_all_replicates_scale_sig), ncol = ncol(Spliceosomal_DEGs_log_z_scaled))
p_beta1_matrix <- matrix(NA, nrow = ncol(PSI_data_all_replicates_scale_sig), ncol = ncol(Spliceosomal_DEGs_log_z_scaled))

# Create a vector storing the green/ white cat data
green_white <- factor(ifelse(grepl("_g_", rownames(PSI_data_all_replicates_scale_sig)), "green", "white"))

# lm(PSI_z_score ~ scale(log(SF_GE)) + cell_type)
for (i in 1:ncol(PSI_data_all_replicates_scale_sig)) {
  for (j in 1:ncol(Spliceosomal_DEGs_log_z_scaled)) {
    # Create a data frame for the linear model
    lm_data <- data.frame(
      PSI = PSI_data_all_replicates_scale_sig[, i], # Column i
      SF_GE = Spliceosomal_DEGs_log_z_scaled[, j], # Column j 
      green_white = green_white
    )
    # lm model
    lm_result <- lm(scale(PSI) ~ scale(SF_GE) + green_white, data = lm_data)
    # Extract the coefficients p-values
    coef_summary <- summary(lm_result)$coefficients
    # beta & p-value
    beta_matrix[i, j] <- coef_summary["scale(SF_GE)", "Estimate"]
    p_beta_matrix[i, j] <- coef_summary["scale(SF_GE)", "Pr(>|t|)"]
    # Store beta1 (coefficient for green_white) and its p-value
    # Note: The actual name in the matrix depends on the factor levels
    green_white_coef_name <- grep("green_white", rownames(coef_summary), value = TRUE)[1]
    if (!is.na(green_white_coef_name)) {
      beta1_matrix[i, j] <- coef_summary[green_white_coef_name, "Estimate"]
      p_beta1_matrix[i, j] <- coef_summary[green_white_coef_name, "Pr(>|t|)"]
    }
  }
}

# Set colnames and rownames for matrices 
# DEGs coefficient matrix
rownames(beta_matrix) <- colnames(PSI_data_all_replicates_scale_si)
colnames(beta_matrix) <- colnames(Spliceosomal_DEGs_log_z_scaled)
# green/white coefficient matrix
rownames(beta1_matrix) <- colnames(PSI_data_all_replicates_scale_si)
colnames(beta1_matrix) <- colnames(Spliceosomal_DEGs_log_z_scaled)
# p value matrix
rownames(p_beta_matrix) <- colnames(PSI_data_all_replicates_scale_si)
colnames(p_beta_matrix) <- colnames(Spliceosomal_DEGs_log_z_scaled)
rownames(p_beta1_matrix) <- colnames(PSI_data_all_replicates_scale_si)
colnames(p_beta1_matrix) <- colnames(Spliceosomal_DEGs_log_z_scaled)

# Calculate FDR matrix for beta 
fdr_beta_matrix <- p_beta_matrix # create fdr vector for SF_GE
for (i in 1:nrow(p_beta_matrix)) {
  fdr_beta_matrix[i,] <- p.adjust(p_beta_matrix[i,], method = "fdr")
}

# Calculate FDR matrix for beta1
fdr_beta1_matrix <- p_beta1_matrix # create fdr vector for cell_type
for (i in 1:nrow(p_beta1_matrix)) {
  fdr_beta1_matrix[i,] <- p.adjust(p_beta1_matrix[i,], method = "fdr")
}

# Filter beta_matrix_rep_1 and beta1_matrix_rep_1 based on FDR 
beta_matrix[fdr_beta_matrix >= 0.1] <- 0
beta1_matrix[fdr_beta1_matrix >= 0.1] <- 0
beta_matrix[beta1_matrix != 0] <- 0

# Remove all 0 rows (introns)
beta_matrix_rep_filt <- beta_matrix_rep_filt[rowSums(is.na(beta_matrix)) != ncol(beta_matrix), ] 

# Calculate mean of all cols 
SF_mean_coef_beta_matrix <- colMeans(beta_matrix, na.rm = TRUE)
SF_mean_coef_beta_matrix <- as.data.frame(SF_mean_coef_beta_matrix)
hist(SF_mean_coef_beta_matrix$SF_mean_coef_beta_matrix)

# SF categories 
SF_mean_coef_beta_matrix <- SF_mean_coef_beta_matrix |>
  mutate(cluster = ifelse(SF_mean_coef_beta_matrix <0, "enhancers", "repressors"))
table(SF_mean_coef_beta_matrix$cluster) # enhancers = 76, repressors = 71

# Create df only contain SF cluster information
col_cluster_SFs <- SF_mean_coef_beta_matrix |>
  rownames_to_column(var = "SF_gene") |>
  dplyr::select(SF_gene, cluster)

# order SF based on mean coefficient 
SF_mean_coef_beta_matrix <- SF_mean_coef_beta_matrix |>
  dplyr::arrange(SF_mean_coef_beta_matrix)

beta_matrix_ordered <- beta_matrix_filt[, rownames(SF_mean_coef_beta_matrix)]
# pheatmap annotation to label splicing enhancers and repressors ? --> create the figure 
anotation_col <- data.frame(Cluster = SF_mean_coef_beta_matrix$cluster)
rownames(anotation_col) <- rownames(SF_mean_coef_beta_matrix)

# Find data limits
max_val <- max(beta_matrix_rep_1_filt, na.rm = TRUE)
min_val <- min(beta_matrix_rep_1_filt, na.rm = TRUE)
abs_limit <- max(abs(min_val), abs(max_val))

# Create symmetric breaks around 0
breaks <- seq(-abs_limit, abs_limit, length.out = 101)  # 100 colors = 101 breaks
View(beta_matrix_ordered)

# Plot
Fig5D <- pheatmap(beta_matrix_ordered,
         show_colnames = FALSE, show_rownames = FALSE, cluster_cols = FALSE,
         annotation_col = anotation_col,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
         breaks = breaks)

## 002. Figure 5E: compare gene expression level of "splicing enhancers" and "splicing repressors" in green and white cells
# Import the avg_normalized_count data 
avg_normalized_count <- read_csv('avg_normalized_count.csv')

# Check the expression of repressors and enhancers in different cell types 
avg_normalized_count_SFs <- avg_normalized_count |>
  filter(gene_id %in% col_cluster_SFs$SF_gene)

# Merge SF cluster information with its GE df 
avg_normalized_count_SFs <- avg_normalized_count_SFs|>
  left_join(col_cluster_SFs, by = c("gene_id" = "SF_gene"))

# pivot longer 
avg_normalized_count_SFs_long <- avg_normalized_count_SFs |>
  pivot_longer(cols = starts_with("avg"), 
               names_to = "sample", 
               values_to = "normalized_counts")
avg_normalized_count_SFs_long_1 <- avg_normalized_count_SFs_long |>
  mutate(cell_type = if_else(grepl("^avg_g_", sample), "green", "white"))

avg_normalized_count_SFs_long_1$cluster <- as.factor(avg_normalized_count_SFs_long_1$cluster)

# Plot 
Fig5E <- ggplot(avg_normalized_count_SFs_long_1, aes(x = cluster, y = log2(normalized_counts + 1), fill = cell_type)) +
  geom_boxplot(width = 0.5, outliers = FALSE) +
  stat_compare_means(aes(group = cell_type), method = "wilcox.test") +
  theme_classic2() +
  scale_fill_manual(values=c("#537D5D", "#999999"))