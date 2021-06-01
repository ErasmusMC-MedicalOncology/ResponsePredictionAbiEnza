# Author:   Job van Riet
# Date:     26-05-21
# Function: Discover biopsy-site associated genes within the full CPCT-02 mCRPC cohort and the Abi/Enza discovery cohort.

set.seed(708)
pacman::p_load('R2CPCT', 'ggplot2', 'DESeq2')

# Required data -----------------------------------------------------------

source('R/3.Figures/misc_Themes.R')

# Load count-data of the full mCRPC dataset (CPCT-02) 
load('/mnt/data2/hartwig/DR71/Oct2020/RData/DR71.RNASeq.RData')

# Load count-data of the discovery cohort.
load('/mnt/data2/hartwig/DR71/Apr2021_AbiEnza/RData/DESeq2Counts.AbiEnza_AllGenes.RData')

# Remove repeated biopsies.
DESeq2Counts.AbiEnza <- DESeq2Counts.AbiEnza[,DESeq2Counts.AbiEnza$responderCategory.DESeq2 != 'Repeat']

# Retrieve batch-effect genes from the full CPCT-02 mCRPC dataset.
diffGenes.Batch <- DR71.RNASeq$DESeq2Results.BetweenMajorBiopsySite %>% 
  dplyr::filter(padj <= 0.05, lfcSE <= 1, abs(log2FoldChange) >= 1) %>% 
  dplyr::distinct(ENSEMBL) %>% 
  dplyr::pull(ENSEMBL)


# Removal of batch effect - Biopsy Site -----------------------------------

plotTSNE <- function(x, colData){
  
  x <- x$Y %>% data.frame()
  x$biopsySite <- colData$biopsySite.Generalized
  
  ggplot2::ggplot(x, ggplot2::aes(x = X1, y = X2, fill = biopsySite)) +
    ggplot2::geom_point(shape = 21, size = 2, color = 'black') +
    ggplot2::labs(x = 't-SNE Dimension 1', y = 't-SNE Dimension 2') +
    ggplot2::scale_fill_manual(name = NULL, values = c('Liver' = '#FF3500', 'Bone' = '#FEFEFE', 'Other' = '#4CA947', 'Lung' = '#9E4CD7', 'Lymph node' = '#0A6C94', 'Soft tissue' = '#EDAEAE'), na.value = 'black', guide = ggplot2::guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
    theme_Job
  
}

## t-SNE of the full cohort ----

countData.VST <- DESeq2::vst(DR71.RNASeq$DESeq2.Raw.allSamples, blind = T)
countData.VST$biopsySite.Generalized <- countData.VST$batchEffectBiopsy

# Perform t-SNE on normalized (VST) counts.
plot.TSNE <- plotTSNE(Rtsne::Rtsne(t(assay(countData.VST)), check_duplicates = F, pca = T, theta = 0.5, perplexity = 15, dims = 2, max_iter = 1E4), colData(countData.VST))
plot.TSNE.Normed <- plotTSNE(Rtsne::Rtsne(t(assay(countData.VST[!rownames(countData.VST) %in% diffGenes.Batch,])), check_duplicates = F, pca = T, theta = 0.5, perplexity = 15, dims = 2, max_iter = 1E4), colData(countData.VST))

## t-SNE of the discovery cohort ----

countData.discovery.VST <- DESeq2::vst(DESeq2Counts.AbiEnza, blind = T)
countData.discovery.VST$biopsySite.Generalized <- ifelse(countData.discovery.VST$biopsySite.Generalized %in% c('Bone', 'Liver', 'Lymph node'), countData.discovery.VST$biopsySite.Generalized, 'Other')

# Perform t-SNE on normalized (VST) counts.
plot.discovery.TSNE <- plotTSNE(Rtsne::Rtsne(t(assay(countData.discovery.VST)), check_duplicates = F, pca = T, theta = 0.5, perplexity = 15, dims = 2, max_iter = 1E4), colData(countData.discovery.VST))
plot.discovery.TSNE.Normed <- plotTSNE(Rtsne::Rtsne(t(assay(countData.discovery.VST[!rownames(countData.discovery.VST) %in% diffGenes.Batch,])), check_duplicates = F, pca = T, theta = 0.5, perplexity = 15, dims = 2, max_iter = 1E4), colData(countData.discovery.VST))

### Combine plots -----------------------------------------------------------

plot.TSNE + plot.TSNE.Normed + 
  plot.discovery.TSNE + plot.discovery.TSNE.Normed +
  patchwork::plot_layout(guides = 'collect', nrow = 2) +
  patchwork::plot_annotation(tag_levels = 'a')
