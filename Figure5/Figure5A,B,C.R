# Load library
library(tidyverse)
library(ggfortify)
library(clusterProfiler)
library(forcats)
library(enrichplot)

## 001. Figure 5A. PCA plot 
# Import PSI data of DSIs for all replicates 
PSI_data_all_replicates_scale_sig <- read_csv("PSI_sigRIs_all_samples_z_scaled.csv")

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

## 002. Figure 5B: DSIs k-means clustering
# Import k-means clustering results 
Sig_RIs_g_w_cluster <- read_csv("Sig_RIs_g_w_cluster.csv")
# Plot
Fig5B <- pheatmap(
  Sig_RIs_g_w_cluster[, 1:16],
  show_rownames = FALSE,
  show_colnames = TRUE,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = data.frame(Cluster = km.res$cluster),
  scale = 'none',
  color = (brewer.pal(10, "RdYlBu")),
  fontsize = 10,
  fontsize_col = 8,
  annotation_colors = list(
    Cluster = colorRampPalette(brewer.pal(8, "Dark2"))(length(unique(km.res$cluster)))
  ),
  annotation_names_row = TRUE,
  annotation_legend = TRUE,
  border_color = NA
)

## 003. Figure 5C: GO terms enrichment analysis
# Import gene ids
intron_id_to_gene_id_sig <- read_csv("intron_id_to_gene_id_sig.csv")
gene_id_introns <- intron_id_to_gene_id_sig$Gene_ID |>
  unique()

# Import GO term data
go_terms_df <- read_csv("go_terms_df.csv")
term2gene <- go_terms_df |>
  dplyr::select(GO, GID)
term2name <- go_terms_df |>
  dplyr::select(GO, TERM)

# GO terms enrichment
BPenrich <- enricher(gene=gene_id_introns,TERM2GENE=term2gene, TERM2NAME = term2name, pAdjustMethod='hochberg',pvalueCutoff=0.05, qvalueCutoff = 0.05)
# mutate enrich ratio = GeneRatio/BgRatio 
y <- mutate(BPenrich, enrich_score = (as.numeric(sub("/.*", "", GeneRatio)) /
                                        as.numeric(sub(".*?/", "", GeneRatio))) /
              (as.numeric(sub("/.*", "", BgRatio)) /
                 as.numeric(sub(".*?/", "", BgRatio))))

# Plot
Fig5C <- ggplot(y, showCategory = 10,
       aes(enrich_score, fct_reorder(Description, enrich_score))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count, shape = Ontology)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() +
  xlab("enrich score") +
  ylab(NULL)

