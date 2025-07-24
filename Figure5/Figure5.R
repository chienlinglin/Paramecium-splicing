# Load library
library(tidyverse)
library(ggfortify)
library(clusterProfiler)
library(forcats)
library(enrichplot)

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

## 002. Figure 5B: DSIs k-means clustering

## 003. Figure 5C: GO terms enrichment analysis
# Import gene ids, term2gene, term2name

term2gene <-
term2name <-

# GO terms enrichment
BPenrich <- enricher(gene=gene_id_introns,TERM2GENE=term2gene, TERM2NAME = term2name, pAdjustMethod='hochberg',pvalueCutoff=0.05, qvalueCutoff = 0.05)
# mutate enrich ratio = GeneRatio/BgRatio 
y <- mutate(BPenrich, enrich_score = (as.numeric(sub("/.*", "", GeneRatio)) /
                                        as.numeric(sub(".*?/", "", GeneRatio))) /
              (as.numeric(sub("/.*", "", BgRatio)) /
                 as.numeric(sub(".*?/", "", BgRatio))))

# Plot
ggplot(y, showCategory = 10,
       aes(enrich_score, fct_reorder(Description, enrich_score))) +
  geom_segment(aes(xend=0, yend = Description)) +
  geom_point(aes(color=p.adjust, size = Count, shape = Ontology)) +
  scale_color_viridis_c(guide=guide_colorbar(reverse=TRUE)) +
  scale_size_continuous(range=c(2, 10)) +
  theme_minimal() +
  xlab("enrich score") +
  ylab(NULL)

