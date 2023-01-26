# Author:    Job van Riet
# Date:      26-01-23
# Function:  Figure of the differential analysis between good vs. poor responders on ARSI-treatment.

# Set seed for reproducibility of t-SNE
set.seed(708813)

# Libraries and data. ----

library(R2CPCT)
library(DESeq2)
library(ggplot2)

# Load ggplot2 themes.
source('R/3.Figures/misc_Themes.R')

# Load metadata of the Abi/Enza-treated patients.
load('/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.Metadata.RData')
load('/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.RNASeq.RData')
DE.LOOCV <- readRDS('/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.LOOCV.Rds')

AbiEnza.Metadata <- AbiEnza.Metadata %>% dplyr::mutate(Responder = ifelse(Responder == 'Bad Responder (≤100 days)', 'Poor Responder (≤100 days)', Responder))

# Function to generate Z-scores.
makeZ <- function (x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m)/s)
}

# List to contain plots.
plots <- list()


# Generate heatmap of DE ----

## Retrieve VST-counts. ----
countData <- SummarizedExperiment::assay(DESeq2::vst(AbiEnza.RNASeq$DESeq2, blind = T))
rownames(countData) <- SummarizedExperiment::rowData(AbiEnza.RNASeq$DESeq2)$SYMBOL

sigGenes <- AbiEnza.RNASeq$DESeq2Results %>% dplyr::filter(AbiEnza.RNASeq$DESeq2Results$isSig == 'Significant')
countData <- countData[rownames(countData) %in% sigGenes$SYMBOL,]

## Cluster rows and columns. ----

orderGenes <- sigGenes %>% dplyr::arrange(log2FoldChange) %>% dplyr::pull(SYMBOL)
orderSamples <- stats::hclust(stats::dist(t(makeZ(countData)), method = 'euclidean'), method = 'ward.D2')

### Add ordering to counts. ----

heatData <- reshape2::melt(makeZ(countData)) %>% 
    dplyr::mutate(
        Var1 = factor(Var1, levels = orderGenes),
        Var2 = factor(Var2, levels = orderSamples$labels[orderSamples$order])
    )

### Generate dendograms. ----

plots$dendro.samples <- ggdendro::ggdendrogram(orderSamples, labels = F) +
    ggplot2::scale_y_discrete(expand=c(0, 0)) +
    ggplot2::scale_x_discrete(expand=c(0.005, 0)) +
    ggplot2::theme(
        axis.text.x = ggplot2::element_blank(),
        axis.text.y = ggplot2::element_blank(),
    )


## Generate annotation tracks. ----

### Treatment duration. ----

plots$TreatmentDuration <- AbiEnza.Metadata %>%
    dplyr::distinct(hmfSampleId, Treatment, treatmentduration_days) %>%
    dplyr::filter(hmfSampleId %in% levels(heatData$Var2)) %>% 
    dplyr::mutate(hmfSampleId = factor(hmfSampleId, levels = levels(heatData$Var2))) %>%
    #Plot.
    ggplot2::ggplot(., ggplot2::aes(x = hmfSampleId, xend = hmfSampleId, yend = 0, y = treatmentduration_days, fill = Treatment)) +
    ggplot2::geom_segment() +
    ggplot2::geom_point(shape = 21, size = 1.25) +
    ggplot2::scale_y_continuous(expand = c(0,0), trans = scales::sqrt_trans(), breaks = c(0, 100, 250, 500, 1000, 2000, 2500), limits = c(0, 2600)) +
    ggplot2::geom_hline(yintercept = 180, color = 'black', lty = 15, lwd = 1) +
    ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(AbiEnza.Metadata$Treatment), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
    ggplot2::labs(y = 'Treatment duration<br><span style = "font-size:5pt">(in days; √)</span>') +
    themeTrack_Job

### Responder class. ----

plots$Responder <- AbiEnza.Metadata %>%
    dplyr::distinct(hmfSampleId, Responder) %>%
    dplyr::filter(hmfSampleId %in% levels(heatData$Var2)) %>% 
    dplyr::mutate(hmfSampleId = factor(hmfSampleId, levels = levels(heatData$Var2))) %>%
    ggplot2::ggplot(., ggplot2::aes(x = hmfSampleId, y = 'Responder class', fill = Responder)) +
    ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
    ggplot2::labs(y = NULL, x = NULL) +
    ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(AbiEnza.Metadata$Responder), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
    themeAnno_Job

### Biopsy sites. ----

plots$Biopsy <- AbiEnza.Metadata %>%
    dplyr::distinct(hmfSampleId, Biopsysite_consolidated) %>%
    dplyr::filter(hmfSampleId %in% levels(heatData$Var2)) %>% 
    dplyr::mutate(hmfSampleId = factor(hmfSampleId, levels = levels(heatData$Var2))) %>%
    #Plot.
    ggplot2::ggplot(., ggplot2::aes(x = hmfSampleId, y = 'Biopsy site', fill = Biopsysite_consolidated)) +
    ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
    ggplot2::labs(y = NULL, x = NULL) +
    ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(AbiEnza.Metadata$Biopsysite_consolidated), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 2, keywidth = 0.5, keyheight = 0.5)) +
    themeAnno_Job


## Log2FC. ----

plots$logFC <- sigGenes %>% 
    dplyr::mutate(SYMBOL = factor(SYMBOL, levels = levels(heatData$Var1))) %>%
    
    ggplot2::ggplot(., aes(x = log2FoldChange, xend = 0, y = SYMBOL, yend = SYMBOL, fill = log2FoldChange)) +
    ggplot2::scale_x_continuous(limits = c(-2.5, 2.5), expand = c(0,0)) +
    
    # Log2FC + Log2SE
    ggplot2::geom_errorbar(aes(xmin = log2FoldChange + lfcSE, xmax = log2FoldChange - lfcSE), width = .75) +
    ggplot2::geom_point(shape = 21) +
    ggplot2::scale_fill_gradient2(limits = c(-2.5, 2.5), low = '#005AB5', mid = 'white', midpoint = 0, high = '#DC3220', guide = ggplot2::guide_colorbar(title = NULL, title.position = 'top', direction = 'horizontal', title.hjust = 0.5, barwidth = 4, barheight = .75)) +
    
    # Guidelines.
    ggplot2::geom_vline(xintercept = 0, lwd = .2, color = 'black', lty = 'dashed') +
    ggplot2::geom_vline(xintercept = .5, lwd = .2, color = 'grey25', lty = 'dotted') +
    ggplot2::geom_vline(xintercept = -.5, lwd = .2, color = 'grey25', lty = 'dotted') +
    
    # Themes.
    ggplot2::labs(x = 'Log<sub>2</sub>FC', y = NULL) +
    theme_Job +
    ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank()
    )


## Heatmap. ----

plots$heatmap <- heatData %>% 
    dplyr::mutate(
        value = ifelse(value > 5, 5, value),
        value = ifelse(value < -5, -5, value)
    ) %>% 
    ggplot2::ggplot(., ggplot2::aes(x = Var2, y = Var1, fill = value))+
    ggplot2::geom_tile(color = 'grey80') +
    ggplot2::labs(x = 'Samples', y = 'DEGs (<i>n</i> = 151)') +
    ggplot2::scale_y_discrete(expand=c(0, 0)) +
    ggplot2::scale_fill_gradient2(limits = c(-5, 5), breaks = c(-5, -2.5, 0, 2.5, 5), labels = c('≤5', -2.5, 0, 2.5, '≥5'), low = '#005AB5', mid = 'white', midpoint = 0, high = '#DC3220', guide = ggplot2::guide_colorbar(title = NULL, title.position = 'top', direction = 'vertical', title.hjust = 0.5, barwidth = .75, barheight = 6)) +
    theme_Job +
    ggplot2::theme(
        text = ggplot2::element_text(size = 7, family='Roboto', face = 'bold'),
        legend.position = 'right',
        axis.text.x = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank()
    )

## Robustness scores. ----

plots$LOOCV <- DE.LOOCV$DE.LOOCV.Summed %>% 
    dplyr::filter(SYMBOL %in% heatData$Var1) %>% 
    dplyr::mutate(SYMBOL = factor(SYMBOL, levels = levels(heatData$Var1))) %>%
    tidyr::complete(SYMBOL, fill = list('totalSig' = 0)) %>% 
    
    ggplot2::ggplot(., aes(x = totalSig, xend = 0, y = SYMBOL, yend = SYMBOL, fill = totalSig)) +
    ggplot2::scale_x_continuous(limits = c(0, 80), expand = c(0,0)) +
    
    ggplot2::geom_segment() +
    ggplot2::geom_point(shape = 21, size = 1, fill = 'black') +
    
    # Themes.
    ggplot2::labs(x = 'LOOCV', y = NULL) +
    theme_Job +
    ggplot2::theme(
        axis.text.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        panel.grid.major.y = ggplot2::element_blank()
    )

# Combine plots. ----

layout <- 'A##
BCD
E##
F##
G##'

svglite::svglite(file = 'Fig4.svg', width = 14, height = 9.5)
plots$TreatmentDuration +
    plots$heatmap + plots$logFC + plots$LOOCV +
    plots$Responder + 
    plots$Biopsy +
    plots$dendro.samples + scale_y_reverse() +
    patchwork::plot_layout(design = layout, guides = 'keep', heights = c(.2, 1, .025, .025, .075), widths = c(1, .2, .11))
dev.off()
