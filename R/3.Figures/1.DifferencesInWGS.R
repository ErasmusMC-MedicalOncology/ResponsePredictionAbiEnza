# Author:    Job van Riet
# Date:      07-01-22
# Function:  Visualize genomic differences.

# Libraries ----

library(R2CPCT)
library(ggplot2)
library(extrafont)
library(patchwork)

# Load themes.
source('R/3.Figures/misc_Themes.R')

# Load metadata of the Abi/Enza-treated patients.
load('/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.Metadata.RData')

# Retrieve WGS-data.
load('/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.Results.RData')

# List to hold plots.
plots <- list()

# Visualization - Mutational burdens ----

## Genome-wide TMB. ----
plots$TMB <- AbiEnza.Results$mutationalBurden %>% 
  dplyr::inner_join(AbiEnza.Metadata, by = c('sample' = 'sampleId')) %>% 
  dplyr::filter(!grepl('Ambi', Responder)) %>% 
  dplyr::group_by(Responder) %>%
  dplyr::mutate(
    median = round(median(Genome.TMB), 2),
    Q1 = round(summary(Genome.TMB)[2], 2),
    Q3 = round(summary(Genome.TMB)[5], 2),
    label = median
    ) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(., ggplot2::aes(x = Responder, y = Genome.TMB, fill = Responder, label = label)) +
  gghalves::geom_half_boxplot(side = 'l', alpha = .3, outlier.shape = NA, notch = F, show.legend = F) +
  gghalves::geom_half_point_panel(side = 'r', position = ggbeeswarm::position_quasirandom(width = .15), size = 1.5, shape = 21, color = 'black', lwd = .01) +
  ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), expand = c(0,0), breaks = c(0:5, 10, 25, 50, 100), limits = c(0, 175)) +
  ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(AbiEnza.Metadata$Responder), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
  ggplot2::stat_summary(fun = min, colour = 'black', geom = 'text', size = 3, show.legend = FALSE, angle = 0) +
  ggplot2::labs(x = NULL, y = 'Tumor Mutational Burden<br><span style = "font-size:5pt">(Genome-wide; log<sub>10</sub>)</span>') +
  ggpubr::geom_bracket(data = AbiEnza.Results$differencesWGS$MutationalBurden %>% dplyr::filter(variable == 'Genome.TMB'), ggplot2::aes(xmin = group1, xmax = group2, label = p.adj.signif, fill = NULL, color = NULL, shape = NULL), y.position = 4.7, tip.length = .01) +
  theme_Job

## Total SV burden. ----
plots$totalSV <- AbiEnza.Results$mutationalBurden %>% 
  dplyr::inner_join(AbiEnza.Metadata, by = c('sample' = 'sampleId')) %>% 
  dplyr::group_by(Responder) %>%
  dplyr::filter(!grepl('Ambi', Responder)) %>% 
  dplyr::mutate(
    median = round(median(totalSV), 0),
    Q1 = round(summary(totalSV)[2], 0),
    Q3 = round(summary(totalSV)[5], 0),
    label = median
  ) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(., ggplot2::aes(x = Responder, y = totalSV, fill = Responder, label = label)) +
  gghalves::geom_half_boxplot(side = 'l', alpha = .3, outlier.shape = NA, notch = F, show.legend = F) +
  gghalves::geom_half_point_panel(side = 'r', position = ggbeeswarm::position_quasirandom(width = .15), size = 1.5, shape = 21, color = 'black', lwd = .01) +
  ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(AbiEnza.Metadata$Responder), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
  ggplot2::stat_summary(fun = min, colour = 'black', geom = 'text', size = 3, show.legend = FALSE, angle = 0) +
  ggplot2::labs(x = NULL, y = 'No. of structural variants<br><span style = "font-size:5pt">(Genome-wide; log10)') +
  ggpubr::geom_bracket(data = AbiEnza.Results$differencesWGS$MutationalBurden %>% dplyr::filter(variable == 'totalSV'), ggplot2::aes(xmin = group1, xmax = group2, label = p.adj.signif, fill = NULL, color = NULL, shape = NULL), y.position = 2000, tip.length = .01) +
  ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), expand = c(0,0), breaks = c(0, 25, 50, 100, 250, 500, 1000, 1500, 2000, 2500), limits = c(25, 2501)) +
  theme_Job

## Total Deletion burden. ----
plots$SV.DEL <- AbiEnza.Results$mutationalBurden %>% 
  dplyr::inner_join(AbiEnza.Metadata, by = c('sample' = 'sampleId')) %>% 
  dplyr::filter(!grepl('Ambi', Responder)) %>% 
  dplyr::group_by(Responder) %>%
  dplyr::mutate(
    median = round(median(SV.DEL), 0),
    Q1 = round(summary(SV.DEL)[2], 0),
    Q3 = round(summary(SV.DEL)[5], 0),
    label = median
  ) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(., ggplot2::aes(x = Responder, y = SV.DEL, fill = Responder, label = label)) +
  gghalves::geom_half_boxplot(side = 'l', alpha = .3, outlier.shape = NA, notch = F, show.legend = F) +
  gghalves::geom_half_point_panel(side = 'r', position = ggbeeswarm::position_quasirandom(width = .15), size = 1.5, shape = 21, color = 'black', lwd = .01) +
  ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), expand = c(0,0), breaks = c(0, 5, 10, 25, 50, 100, 250, 500, 1000, 1500), limits = c(0, 1501)) +
  ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(AbiEnza.Metadata$Responder), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
  ggplot2::stat_summary(fun = min, colour = 'black', geom = 'text', size = 3, show.legend = FALSE, angle = 0) +
  ggplot2::labs(x = NULL, y = 'No. of Deletions<br><span style = "font-size:5pt">(Genome-wide; log10)') +
  ggplot2::scale_shape_manual(values = c(21, 23)) +
  ggpubr::geom_bracket(data = AbiEnza.Results$differencesWGS$MutationalBurden %>% dplyr::filter(variable == 'SV.DEL'), ggplot2::aes(xmin = group1, xmax = group2, label = p.adj.signif, fill = NULL, color = NULL, shape = NULL), y.position = 1000, tip.length = .01) +
  theme_Job

## Total Tandem duplication burden. ----
plots$SV.DUP <- AbiEnza.Results$mutationalBurden %>% 
  dplyr::inner_join(AbiEnza.Metadata, by = c('sample' = 'sampleId')) %>% 
  dplyr::filter(!grepl('Ambi', Responder)) %>% 
  dplyr::group_by(Responder) %>%
  dplyr::mutate(
    median = round(median(SV.DUP), 0),
    Q1 = round(summary(SV.DUP)[2], 0),
    Q3 = round(summary(SV.DUP)[5], 0),
    label = median
  ) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(., ggplot2::aes(x = Responder, y = SV.DUP, fill = Responder, label = label)) +
  gghalves::geom_half_boxplot(side = 'l', alpha = .3, outlier.shape = NA, notch = F, show.legend = F) +
  gghalves::geom_half_point_panel(side = 'r', position = ggbeeswarm::position_quasirandom(width = .15), size = 1.5, shape = 21, color = 'black', lwd = .01) +
  ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), expand = c(0,0), breaks = c(0, 5, 10, 25, 50, 100, 250, 500, 1000, 1500), limits = c(0, 1501)) +
  ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(AbiEnza.Metadata$Responder), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
  ggplot2::stat_summary(fun = min, colour = 'black', geom = 'text', size = 3, show.legend = FALSE, angle = 0) +
  ggplot2::labs(x = NULL, y = 'No. of Tandem Duplications<br><span style = "font-size:5pt">(Genome-wide; log10)') +
  ggplot2::scale_shape_manual(values = c(21, 23)) +
  ggpubr::geom_bracket(data = AbiEnza.Results$differencesWGS$MutationalBurden %>% dplyr::filter(variable == 'SV.DUP'), ggplot2::aes(xmin = group1, xmax = group2, label = p.adj.signif, fill = NULL, color = NULL, shape = NULL), y.position = 1000, tip.length = .01) +
  theme_Job


### Combine plots. ----

plots$TMB + plots$totalSV + plots$SV.DUP + plots$SV.DEL +
  patchwork::plot_layout(guides = 'collect', nrow = 1) +
  patchwork::plot_annotation(tag_levels = 'a')


# Visualization - Chromosomal Arms ----

data.Arms <- AbiEnza.Results$GISTIC2.AllSamples$gisticBroadScoresPerArm %>% reshape2::dcast(variable ~ `Chromosome Arm`)
rownames(data.Arms) <- data.Arms$variable; data.Arms$variable <- NULL

data.Arms <- data.Arms[rownames(data.Arms) %in% (AbiEnza.Metadata %>% dplyr::filter(!grepl('Ambi', Responder)) %>% dplyr::pull(sampleId)),]

# Column annotation.
annotation.col <- data.frame(
  'Responder class' = AbiEnza.Metadata$Responder,
  'Biopsy site (Generalized)' = AbiEnza.Metadata$Biopsysite_consolidated,
  'Treatment (Generalized)' = AbiEnza.Metadata$Treatment,
  row.names = AbiEnza.Metadata$sampleId, check.names = F
)

# Colors of the annotations.
annotation.colors <- list(
  'Responder class' = colorPalette[names(colorPalette) %in% unique(annotation.col$`Responder class`)],
  'Biopsy site (Generalized)' = colorPalette[names(colorPalette) %in% unique(annotation.col$`Biopsy site (Generalized)`)],
  'Treatment (Generalized)' = colorPalette[names(colorPalette) %in% unique(annotation.col$`Treatment (Generalized)`)]
)

pheatmap::pheatmap(
  t(data.Arms), 
  fontsize_row = 6,  treeheight_col = 15, treeheight_row = 15,
  angle_col = 90,
  show_colnames = F, show_rownames = T,
  cluster_cols = T, cluster_rows = F, scale = 'row', fontsize = 8,
  annotation_col = annotation.col, annotation_colors = annotation.colors,
  clustering_distance_rows = 'euclidean', clustering_distance_cols = 'euclidean', clustering_method = 'ward.D2',
  cellheight = 6, cellwidth = 3,
  border_color = 'grey90', drop_levels = T, 
  color = colorRampPalette(c('darkgreen', 'green', 'white', 'red', 'darkred'))(101)
)
