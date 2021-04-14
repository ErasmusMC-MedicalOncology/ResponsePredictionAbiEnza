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
AbiEnza.Metadata <- readxl::read_xlsx('Misc/Suppl. Table 1 - OverviewOfData.xlsx', sheet = 'Sample overview')

# Load DE-Genes from DESeq2.
AbiEnza.DE <- readxl::read_xlsx('Misc/Suppl. Table 1 - OverviewOfData.xlsx', sheet = 'Differential Expression') %>% dplyr::filter(`Significant Threshold` == 'Significant')

# Retrieve RNA-Seq counts.
load('/mnt/data2/hartwig/DR71/Oct2020_AbiEnza/RData/DESeq2Counts.AbiEnza.RData')

# List to contain plots.
plots <- list()


# Generate t-SNE of Poor vs. Good responders ------------------------------

## Retrieve VST-counts. ----
countData <- SummarizedExperiment::assay(DESeq2::vst(DESeq2Counts.AbiEnza, blind = F))
rownames(countData) <- rowData(DESeq2Counts.AbiEnza)$SYMBOL

## Perform T-SNE on DE-genes. ----
TSNE.AbiEnza.All <- Rtsne::Rtsne(t(countData), check_duplicates = T, pca = T, theta = .5, perplexity = 9, dims = 2, max_iter = 1E5, num_threads = 20)
TSNE.AbiEnza.All <- TSNE.AbiEnza.All$Y %>% data.frame()
TSNE.AbiEnza.All$responderCategory <- DESeq2Counts.AbiEnza$responderCategory
TSNE.AbiEnza.All$ClassificationType <- DESeq2Counts.AbiEnza$ClassificationType
TSNE.AbiEnza.All$Sample <- DESeq2Counts.AbiEnza$sampleId

## Perform T-SNE on DE-genes. ----
TSNE.AbiEnza <- Rtsne::Rtsne(t(countData[rownames(countData) %in% AbiEnza.DE$SYMBOL,]), check_duplicates = T, pca = T, theta = .5, perplexity = 9, dims = 2, max_iter = 1E5, num_threads = 20)
TSNE.AbiEnza <- TSNE.AbiEnza$Y %>% data.frame()
TSNE.AbiEnza$responderCategory <- DESeq2Counts.AbiEnza$responderCategory
TSNE.AbiEnza$ClassificationType <- DESeq2Counts.AbiEnza$ClassificationType
TSNE.AbiEnza$Sample <- DESeq2Counts.AbiEnza$sampleId

## Plot t-SNE ----

plots$tSNE.All <- ggplot2::ggplot(TSNE.AbiEnza.All, ggplot2::aes(x = X1, y = X2, fill = responderCategory, shape = ClassificationType, label = Sample)) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::scale_fill_manual(values = c('Poor Responder (≤100 days)' = '#E69F00', 'Good Responder (>100 days)' = '#019E73', 'Unknown Responder' = '#999999')) +
  ggplot2::scale_shape_manual(values = c(21, 23)) +
  ggplot2::labs(x = 't-SNE Dimension 1', y = 't-SNE Dimension 2') +
  ggplot2::guides(fill = ggplot2::guide_legend(title = 'Responder Category', title.position = 'top', title.hjust = 0.5, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
  theme_Job

plots$tSNE.DE <- ggplot2::ggplot(TSNE.AbiEnza, ggplot2::aes(x = X1, y = X2, fill = responderCategory, shape = ClassificationType, label = Sample)) +
  ggplot2::geom_point(size = 2.5) +
  ggplot2::scale_fill_manual(values = c('Poor Responder (≤100 days)' = '#E69F00', 'Good Responder (>100 days)' = '#019E73', 'Unknown Responder' = '#999999')) +
  ggplot2::scale_shape_manual(values = c(21, 23)) +
  ggplot2::labs(x = 't-SNE Dimension 1', y = 't-SNE Dimension 2') +
  ggplot2::guides(fill = ggplot2::guide_legend(title = 'Responder Category', title.position = 'top', title.hjust = 0.5, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
  theme_Job


# Generate heatmap --------------------------------------------------------

# Heatmap.
heatData <- countData[rownames(countData) %in% AbiEnza.DE$SYMBOL,]

# Column annotation.
annotation.row <- data.frame('Direction' = factor(ifelse(colnames(t(heatData)) %in% (AbiEnza.DE %>% dplyr::filter(`log2FoldChange (Poor vs. Good)` >= 0))$SYMBOL, 'Up-regulated in Poor Responders', 'Down-regulated in Poor Responders')), row.names = colnames(t(heatData)))
annotation.col <- data.frame(
  'Responder.Category' = DESeq2Counts.AbiEnza$responderCategory,
  'TMPRSS2.ERG' = DESeq2Counts.AbiEnza$hasGenomicERG,
  'Subtype' = DESeq2Counts.AbiEnza$ClassificationType,
  'Treatment' = DESeq2Counts.AbiEnza$treatment,
  'Treatment.Duration' = DESeq2Counts.AbiEnza$treatmentDurationInDays,
  'Biopsy.Site' = DESeq2Counts.AbiEnza$biopsySite.Generalized,
  row.names = DESeq2Counts.AbiEnza$sampleId
)

# Clean-up annotations.
annotation.col <- annotation.col %>% dplyr::mutate(Biopsy.Site = ifelse(Biopsy.Site %in% c('Bone', 'Lung', 'Lymph node', 'Liver', 'Prostate'), Biopsy.Site, 'Other'))
annotation.col <- annotation.col %>% dplyr::mutate(Treatment = gsub('/.*', '', Treatment))
annotation.col <- annotation.col %>% dplyr::mutate(Treatment.Duration = ifelse(is.na(Treatment.Duration), -50, Treatment.Duration))

# Colors of the annotations.
annotation.colors <- list(
  'Direction' = c('Up-regulated in Poor Responders' = '#D03C3F', 'Down-regulated in Poor Responders' = '#5EA153'),
  'Responder.Category' = c('Poor Responder (≤100 days)' = '#E69F00', 'Good Responder (>100 days)' = '#019E73', 'Unknown Responder' = '#999999'),
  'Biopsy.Site' = c('Liver' = '#FF3500', 'Lung' = '#FFA000', 'Prostate' = '#EDAEAE', 'Bone' = '#FEFEFE', 'Other' = '#4CA947', 'Lymph node' = '#0A6C94'),
  'TMPRSS2.ERG' = c('Yes' = 'grey10', 'No' = 'grey80'),
  'Treatment' = c('Abiraterone' = '#2a7fff', 'Enzalutamide' = '#ff7f2a'),
  'Subtype' = c('t-NEPC (RNA-Seq)' = 'grey10', 'mCRPC' = 'grey80'),
  'Treatment.Duration' = RColorBrewer::brewer.pal(9, 'Greens')
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
  patchwork::plot_layout(guides = 'keep', nrow = 1) +
  patchwork::plot_annotation(tag_levels = 'a')
