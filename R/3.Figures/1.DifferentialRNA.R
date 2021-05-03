# Author:    Job van Riet
# Date:      08-04-21
# Function:  Figure of the differential analysis between good vs. poor responders on Abi/Enza-treatment.

# Set seed for reproducibility of t-SNE
set.seed(708813)

# Libraries ---------------------------------------------------------------

library(R2CPCT)
library(DESeq2)

# Load ggplot2 themes.
source('R/3.Figures/misc_Themes.R')

# Load metadata of the Abi/Enza-treated patients.
load('/mnt/data2/hartwig/DR71/Apr2021_AbiEnza/RData/AbiEnza.Metadata.RData')

# Load DE-Genes from DESeq2.
AbiEnza.DE <- readxl::read_xlsx('Misc/Suppl. Table 1 - OverviewOfData.xlsx', sheet = 'Differential Expression') %>% dplyr::filter(`Significant Threshold` == 'Significant')

# Retrieve RNA-Seq counts.
load('/mnt/data2/hartwig/DR71/Apr2021_AbiEnza/RData/DESeq2Counts.AbiEnza.RData')

# List to contain plots.
plots <- list()


# Generate t-SNE of Poor vs. Good responders ------------------------------

## Retrieve VST-counts. ----
countData <- SummarizedExperiment::assay(DESeq2::vst(DESeq2Counts.AbiEnza, blind = F))
rownames(countData) <- rowData(DESeq2Counts.AbiEnza)$SYMBOL

## Perform T-SNE on DE-genes. ----
TSNE.AbiEnza.All <- Rtsne::Rtsne(t(countData), check_duplicates = T, pca = T, theta = .5, perplexity = 9, dims = 2, max_iter = 1E5, num_threads = 20)
TSNE.AbiEnza.All <- TSNE.AbiEnza.All$Y %>% data.frame()
TSNE.AbiEnza.All$Responder <- DESeq2Counts.AbiEnza$Responder
TSNE.AbiEnza.All$Sample <- DESeq2Counts.AbiEnza$sampleId

## Perform T-SNE on DE-genes. ----
TSNE.AbiEnza <- Rtsne::Rtsne(t(countData[rownames(countData) %in% AbiEnza.DE$SYMBOL,]), check_duplicates = T, pca = T, theta = .5, perplexity = 9, dims = 2, max_iter = 1E5, num_threads = 20)
TSNE.AbiEnza <- TSNE.AbiEnza$Y %>% data.frame()
TSNE.AbiEnza$Responder <- DESeq2Counts.AbiEnza$Responder
TSNE.AbiEnza$Sample <- DESeq2Counts.AbiEnza$sampleId

## Plot t-SNE ----

plots$tSNE.All <- ggplot2::ggplot(TSNE.AbiEnza.All, ggplot2::aes(x = X1, y = X2, fill = Responder, label = Sample)) +
  ggplot2::geom_point(shape = 21, size = 2.5) +
  ggplot2::scale_fill_manual(values = c('Bad Responder (≤100 days)' = '#E69F00', 'Good Responder (>100 days)' = '#019E73')) +
  ggplot2::scale_shape_manual(values = c(21, 23)) +
  ggplot2::labs(x = 't-SNE Dimension 1', y = 't-SNE Dimension 2') +
  ggplot2::guides(fill = ggplot2::guide_legend(title = 'Responder Category', title.position = 'top', title.hjust = 0.5, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
  theme_Job

plots$tSNE.DE <- ggplot2::ggplot(TSNE.AbiEnza, ggplot2::aes(x = X1, y = X2, fill = Responder, label = Sample)) +
  ggplot2::geom_point(shape = 21, size = 2.5) +
  ggplot2::scale_fill_manual(values = c('Bad Responder (≤100 days)' = '#E69F00', 'Good Responder (>100 days)' = '#019E73')) +
  ggplot2::scale_shape_manual(values = c(21, 23)) +
  ggplot2::labs(x = 't-SNE Dimension 1', y = 't-SNE Dimension 2') +
  ggplot2::guides(fill = ggplot2::guide_legend(title = 'Responder Category', title.position = 'top', title.hjust = 0.5, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
  theme_Job


# Generate heatmap --------------------------------------------------------

# Heatmap.
heatData <- countData[rownames(countData) %in% AbiEnza.DE$SYMBOL,]

# Column annotation.
annotation.row <- data.frame('Direction' = factor(ifelse(colnames(t(heatData)) %in% (AbiEnza.DE %>% dplyr::filter(`log2FoldChange (Poor vs. Good)` >= 0))$SYMBOL, 'Up-regulated in Bad Responders', 'Down-regulated in Bad Responders')), row.names = colnames(t(heatData)))
annotation.col <- data.frame(
  'Responder.category' = DESeq2Counts.AbiEnza$Responder,
  'ERG.fusion' = DESeq2Counts.AbiEnza$hasGenomicERG,
  'Treatment' = DESeq2Counts.AbiEnza$treatment.Generalized,
  'Treatment.duration' = DESeq2Counts.AbiEnza$treatmentDurationInDays,
  'Biopsy.site' = DESeq2Counts.AbiEnza$biopsySite.Generalized,
  row.names = DESeq2Counts.AbiEnza$sampleId
)

# Clean-up annotations.
annotation.col <- annotation.col %>% dplyr::mutate(Treatment = gsub('/.*', '', Treatment))
annotation.col <- annotation.col %>% dplyr::mutate(Treatment.duration = ifelse(is.na(Treatment.duration), -50, Treatment.duration))

# Colors of the annotations.
annotation.colors <- list(
  'Direction' = c('Up-regulated in Bad Responders' = '#D03C3F', 'Down-regulated in Bad Responders' = '#5EA153'),
  'Responder.category' = c('Bad Responder (≤100 days)' = '#E69F00', 'Good Responder (>100 days)' = '#019E73'),
  'Biopsy.site' = c('Liver' = '#FF3500', 'Bone' = '#FEFEFE', 'Other' = '#4CA947', 'Lung' = '#9E4CD7', 'Lymph node' = '#0A6C94', 'Soft tissue' = '#EDAEAE'),
  'ERG.fusion' = c('Yes' = '#2F385E', 'No' = '#FFFFFF'),
  'Treatment' = c('Abiraterone' = '#2a7fff', 'Enzalutamide' = '#ff7f2a'),
  'Treatment.duration' = RColorBrewer::brewer.pal(9, 'Greens')
)

# Plot heatmap.
pheatmap::pheatmap(
  heatData,
  fontsize_row = 6, treeheight_col = 15, treeheight_row = 15,
  angle_col = 90,
  show_colnames = F, show_rownames = T,
  cluster_cols = T, cluster_rows = T, scale = 'row',
  annotation_col = annotation.col, annotation_row = annotation.row, annotation_colors = annotation.colors,
  clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'ward.D2',
  cellheight = 6, cellwidth = 5,
  border_color = 'grey90',
  color = colorRampPalette(c('green', 'white', 'red'))(101), 
  cutree_rows = 2, cutree_cols = 2
)


# Combine plots -----------------------------------------------------------

plots$tSNE.All + plots$tSNE.DE + 
  patchwork::plot_layout(guides = 'collect', nrow = 2) +
  patchwork::plot_annotation(tag_levels = 'a')
