# Author:    Job van Riet
# Date:      21-03-22
# Function:  Detection of AR-V7.

# Libraries ----

library(R2CPCT)
library(patchwork)

# Load metadata of the Abi/Enza-treated patients.
load('/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.Metadata.RData')
load('/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.RNASeq.RData')

# Retrieve WTS of DR-071.
DR71.ARvs <- readRDS('/mnt/share1/repository/HMF/DR71/Dec2021/RData/DR71.ARvs.Rds')

# Load ggplot2 themes.
source('R/3.Figures/misc_Themes.R')

# GSEA ----

plot.GSEA <- AbiEnza.RNASeq$GSEA %>% 
    dplyr::filter(padj <= 0.05) %>% 
    dplyr::arrange(NES) %>% 
    dplyr::inner_join(read.delim('Misc/pathwayClean.csv', sep = '\t'), by = c('pathway' = 'pathwayOriginal')) %>% 
    dplyr::distinct(pathwayClean, NES) %>% 
    dplyr::mutate(
        pathwayClean = gsub(' \\(Hall.*', ' <sup style="color:#4987BA">(H)</sup>', pathwayClean),
        pathwayClean = gsub(' \\(Wiki.*', ' <sup style="color:#F17F64">(W)</sup>', pathwayClean)
    ) %>% 
    ggplot2::ggplot(aes(x = reorder(pathwayClean, NES), xend = reorder(pathwayClean, NES), yend = 0, y = NES, fill = NES)) +
    ggplot2::geom_segment() +
    ggplot2::geom_point(shape = 21, size = 2) +
    ggplot2::geom_hline(yintercept = 0, lty = '11') +
    ggplot2::scale_y_continuous(limits = c(-3, 3)) +
    ggplot2::scale_fill_gradient2(limits = c(-3, 3), breaks = c(-3, 0, 3), low = '#005AB5', mid = 'white', midpoint = 0, high = '#DC3220', guide = ggplot2::guide_colorbar(title = NULL, title.position = 'top', direction = 'horizontal', title.hjust = 0.5, barwidth = 6, barheight = .5)) +
    ggplot2::scale_alpha_continuous(guide = 'none') +
    ggplot2::labs(x = 'GSEA<br>(Hallmark & WikiPathways; <i>q</i> ≤ 0.05)', y = 'Normalized Enrichment Score<br>Poor vs. Good responders') + 
    theme_Job + 
    ggplot2::theme(
        axis.text.y = ggtext::element_markdown(size = 8)
    ) + coord_flip()


# Visualize and test AR-V7 expression. ----

AR.PSI <- DR71.ARvs$AR.PSI %>% 
    dplyr::select(sample, PSI.ARv7, counts.ARv7) %>% 
    dplyr::inner_join(AbiEnza.Metadata, by = c('sample' = 'sampleId')) %>% 
    dplyr::mutate(
        Responder = ifelse(Responder == 'Bad Responder (≤100 days)', 'Poor Responder (≤100 days)', Responder),
        PSI.ARv7 = ifelse(is.nan(PSI.ARv7), 0, PSI.ARv7)
    ) %>% 
    dplyr::group_by(Responder) %>% 
    dplyr::mutate(median = median(PSI.ARv7)) %>% 
    dplyr::ungroup()

# Test.
stat.test <- rstatix::pairwise_wilcox_test(AR.PSI, PSI.ARv7 ~ Responder)

# Plot.
plot.ARv <- ggplot2::ggplot(AR.PSI, ggplot2::aes(x = reorder(Responder, -median), y = PSI.ARv7, fill = Responder)) +
    gghalves::geom_half_boxplot(outlier.shape = NA) +
    gghalves::geom_half_point_panel(shape = 21) +
    ggplot2::scale_y_continuous(labels = scales::percent_format(accuracy = 1), breaks = c(seq(0, .1, .01)), expand = c(0,.001), limits = c(0, .08)) +
    ggplot2::labs(x = NULL, y = 'Rel. expression of variant-specific SJ<br>(SJ<sub><i>x</i></sub> / SJ<sub>Exon1_Exon2</sub>)') +
    ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(AR.PSI$Responder), guide = ggplot2::guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
    ggpubr::geom_bracket(data = stat.test, ggplot2::aes(xmin = group1, xmax = group2, label = p.adj), inherit.aes = F, label.size = 2, y.position = 0.07, step.increase = .05, tip.length = .01) +
    theme_Job +
    ggplot2::theme(panel.grid.minor.y = ggplot2::element_blank())


svglite::svglite(file = 'SupplFig2.svg', width = 9, height = 5)
    plot.GSEA + plot.ARv + patchwork::plot_layout(widths = c(1, 1.25)) + patchwork::plot_annotation(tag_levels = 'a')
dev.off()
