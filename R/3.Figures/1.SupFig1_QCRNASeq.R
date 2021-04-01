# Author:    Job van Riet
# Date:      01-04-21
# Function:  Import and processing of the Abi/Enza-treated RNA-Seq samples (DR-071).

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


# Generate t-SNE of Poor vs. Good responders ------------------------------

# Retrieve VST-counts.
countData <- SummarizedExperiment::assay(DESeq2::vst(DESeq2Counts.AbiEnza, blind = T))
rownames(countData) <- rowData(DESeq2Counts.AbiEnza)$SYMBOL

# Perform T-SNE on all genes.
TSNE.AbiEnza <- Rtsne::Rtsne(t(countData[rownames(countData) %in% AbiEnza.DE$SYMBOL,]), check_duplicates = F, pca = T, theta = 0.5, perplexity = 4, dims = 2, max_iter = 1E4)
TSNE.AbiEnza <- TSNE.AbiEnza$Y %>% data.frame()
TSNE.AbiEnza$Classification <- DESeq2Counts.AbiEnza$responderCategory
TSNE.AbiEnza$Sample <- DESeq2Counts.AbiEnza$sampleId

# Plot the NE-TSNE
ggplot2::ggplot(TSNE.AbiEnza, aes(x = X1, y = X2, fill = Classification, label = Sample)) +
  ggplot2::geom_point(shape = 21, size = 2.5) +
  ggplot2::scale_fill_manual(values = c('Poor Responder (≤100 days)' = '#ef233c', 'Good Responder (>100 days)' = '#55a630', 'Unknown Responder' = '#ffc145')) +
  ggplot2::labs(x = 't-SNE Dimension 1', y = 't-SNE Dimension 2') +
  ggplot2::guides( fill = guide_legend(title = 'Mutational Categories', title.position = 'top', title.hjust = 0.5, ncol = 3, keywidth = 0.5, keyheight = 0.5)) +
  theme_Job
