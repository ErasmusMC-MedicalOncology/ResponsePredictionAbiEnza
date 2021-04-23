# Author:    Job van Riet
# Date:      16-04-21
# Function:  Genomic overview of the Abi/Enza-treated patients.


# Import data and themes --------------------------------------------------

# Load ggplot2 themes.
source('R/3.Figures/misc_Themes.R')

# Load metadata of the Abi/Enza-treated patients.
AbiEnza.Metadata <- readxl::read_xlsx('Misc/Suppl. Table 1 - OverviewOfData.xlsx', sheet = 'Sample overview')

# Load WGS data.
load('/mnt/data2/hartwig/DR71/Oct2020_AbiEnza/RData/AbiEnza.Results.RData')


# Main Figure - Genomic Landscape -----------------------------------------

tracks.landscape <- list()

# Sort on Abi/Enza treatment duration.
orderSamples <- AbiEnza.Metadata %>%
  dplyr::arrange(responderCategory, -treatmentDurationInDays) %>%
  dplyr::pull(sampleId)


## Track - Responder group ----

tracks.landscape$responderCategory <- AbiEnza.Metadata %>%
  dplyr::distinct(sampleId, responderCategory) %>%
  dplyr::mutate(sampleId = factor(sampleId, levels = orderSamples)) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = 'Responder category', fill = responderCategory)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = c('Poor Responder (≤100 days)' = '#E69F00', 'Good Responder (>100 days)' = '#019E73', 'Unknown Responder' = '#999999'), na.value = 'white', guide = ggplot2::guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
  themeAnno_Job


## Track - Treatment duration ----

tracks.landscape$treatmentDuration <- AbiEnza.Metadata %>%
  dplyr::distinct(sampleId, treatment, treatmentDurationInDays) %>%
  dplyr::mutate(sampleId = factor(sampleId, levels = orderSamples)) %>%
  dplyr::mutate(treatment = gsub('/.*', '', treatment)) %>%
  #Plot.
  ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = treatmentDurationInDays, fill = treatment)) +
  ggplot2::geom_bar(stat = 'identity', color = 'black', lwd = .33, width = .8) +
  ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100, 200, 300, 400, 500, 600), limits = c(0, 601)) +
  ggplot2::geom_hline(yintercept = 100, color = 'red', lty = 'dotted', lwd = .5) +
  ggplot2::scale_fill_manual('Treatment', values = c('Abiraterone' = '#2a7fff', 'Enzalutamide' = '#ff7f2a'), guide = ggplot2::guide_legend(title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  ggplot2::labs(y = 'Treatment and duration<br><span style = "font-size:5pt">(in nr. of days)</span>') +
  themeTrack_Job


## Track - Genome-wide TMB ----

tracks.landscape$TMB <- AbiEnza.Results$mutationalBurden %>%
  dplyr::distinct(sample, Genome.TMB) %>%
  dplyr::mutate(sample = factor(sample, levels = orderSamples)) %>%
  #Plot.
  ggplot2::ggplot(., ggplot2::aes(x = sample, y = Genome.TMB)) +
  ggplot2::geom_bar(stat = 'identity', fill = '#F75658', color = 'black', lwd = .33, width = .8) +
  ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(0, 2.5, 5, 10, 25, 50, 100), limits = c(0, 115)) +
  ggplot2::geom_hline(yintercept = 10, color = 'black', lty = 'dotted', lwd = .5) +
  ggplot2::labs(y = 'Tumor Mutational Burden<br><span style = "font-size:5pt">(Genome-wide; log<sub>10</sub>)</span>') +
  ggplot2::coord_trans(y = scales::pseudo_log_trans()) +
  themeTrack_Job

## Track - Overview of structural variants. ----

tracks.landscape$SV <-  AbiEnza.Results$mutationalBurden %>%
  dplyr::distinct(sample, totalSV) %>%
  dplyr::mutate(sample = factor(sample, levels = orderSamples)) %>%
  #Plot.
  ggplot2::ggplot(., ggplot2::aes(x = sample, y = totalSV)) +
  ggplot2::geom_bar(stat = 'identity', fill = '#00AF66', color = 'black', lwd = .33, width = .8) +
  ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(0, 100, 250, 500, 1000, 1500, 2000), limits = c(0, 2100)) +
  ggplot2::labs(y = 'Nr. of structural variants<br><span style = "font-size:5pt">(Genome-wide)') +
  themeTrack_Job


## Track - Overview of structural variant categories. ----

tracks.landscape$SV.Cat <- AbiEnza.Results$mutationalBurden %>%
  dplyr::select(sample, totalSV, SV.TRA, SV.DEL, SV.DUP, SV.SINGLE, SV.INV, SV.INS) %>%
  reshape2::melt(id.vars = c('sample', 'totalSV')) %>%
  dplyr::mutate(
    type = 'Other',
    type = ifelse(variable == 'SV.DEL', 'Deletions', type),
    type = ifelse(variable == 'SV.INS', 'Insertions', type),
    type = ifelse(variable == 'SV.INV', 'Inversions', type),
    type = ifelse(variable == 'SV.DUP', 'Tandem Duplications', type),
    type = ifelse(variable == 'SV.TRA', 'Translocations', type),
    type = ifelse(variable == 'SV.SINGLE', 'Singles', type),
    type = factor(type, levels = c('Deletions', 'Inversions', 'Tandem Duplications', 'Translocations', 'Insertions', 'Singles')),
    value.Rel = value / totalSV,
    value.Rel = ifelse(is.na(value.Rel), 0, value.Rel),
    sample = factor(sample, levels = orderSamples)
  ) %>%
  #Plot.
  ggplot2::ggplot(., ggplot2::aes(x = sample, y = value.Rel, fill = type)) +
  ggplot2::geom_bar(stat = 'identity', color = 'black', lwd = .33, width = .8) +
  ggplot2::scale_y_continuous(expand = c(0,0), labels = scales::percent, breaks = c(0, .25, .5, .75, 1)) +
  ggplot2::scale_fill_manual('Structural variant categories', values = c('Translocations' = '#375D96', 'Deletions' = '#ff8c00', 'Tandem Duplications' ='#fc6769', 'Insertions' = 'yellow', 'Inversions' = 'skyblue', 'Singles' = 'grey75'), guide = ggplot2::guide_legend(title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  ggplot2::labs(y = 'Structural variants<br><span style = "font-size:5pt">Relative categories</span>') +
  themeTrack_Job

## Track - Genome-wide ploidy. ----

tracks.landscape$Ploidy <- AbiEnza.Results$mutationalBurden %>%
  dplyr::select(sample, genomePloidy, wholeGenomeDuplication) %>%
  dplyr::mutate(
    sample = factor(sample, levels = orderSamples),
    wholeGenomeDuplication = ifelse(wholeGenomeDuplication, '*', '')
  ) %>%
  #Plot.
  ggplot2::ggplot(., ggplot2::aes(x = sample, y = genomePloidy, label = wholeGenomeDuplication, fill = genomePloidy)) +
  ggplot2::geom_bar(stat = 'identity', lwd = .33, color = 'black', width = .8) +
  ggplot2::scale_fill_gradient2('Genome-wide ploidy', limits = c(0, 6.5), low = 'darkred', mid = 'white', midpoint = 2, high = 'darkgreen', guide = ggplot2::guide_colorbar(title = NULL, title.position = 'top', direction = 'vertical', title.hjust = 0.5, barwidth = .75, barheight = 3)) +
  ggplot2::scale_y_continuous(expand = c(0,0), limits = c(0, 7)) +
  ggplot2::geom_text() +
  ggplot2::labs(y = 'Genome Ploidy<br><span style = "font-size:5pt">') +
  themeTrack_Job


## Track - Mutational Signatures ----

tracks.landscape$MutSigs.SNV <- R2CPCT::plotMutationalSignaturesCOSMIC(mutSigs = AbiEnza.Results$mutSigs$SNV, orderSamples =  as.character(orderSamples), minContrib = 5, combineSigs = T) + themeTrack_Job


## Track - MSI ----

tracks.landscape$MSI <- AbiEnza.Results$mutationalBurden %>%
  dplyr::select(sample, msStatus) %>%
  dplyr::mutate(
    sample = factor(sample, levels = orderSamples)
  ) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sample, y = 'MSI', fill = msStatus)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = c('MSI' = '#B66A29', 'MSS' = 'white')) +
  themeAnno_Job + ggplot2::theme(legend.position = 'none')


## Track - CHORD. ----

tracks.landscape$HRD <- AbiEnza.Results$mutationalBurden %>%
  dplyr::select(sample, hr_status) %>%
  dplyr::mutate(
    sample = factor(sample, levels = orderSamples),
    isHRD = ifelse(hr_status == 'HR_deficient', 'HR-deficient', 'HR-proficient')
  ) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sample, y = 'HRD', fill = isHRD)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = c('HR-deficient' = '#FF0080', 'HR-proficient' = 'white')) +
  themeAnno_Job + ggplot2::theme(legend.position = 'none')


## Track - Chromothripsis. ----

tracks.landscape$Chromothripsis <- AbiEnza.Results$mutationalBurden %>%
  dplyr::select(sample, hasChromothripsis) %>%
  dplyr::mutate(
    sample = factor(sample, levels = orderSamples)
  ) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sample, y = 'Chromothripsis', fill = hasChromothripsis)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = c('Yes' = '#E64F7050', 'No' = 'white')) +
  themeAnno_Job + ggplot2::theme(legend.position = 'none')


## Track - Biopsy site. ----

tracks.landscape$biopsySite <- AbiEnza.Metadata %>%
  dplyr::select(sampleId, biopsySite.Generalized) %>%
  dplyr::mutate(
    sampleId = factor(sampleId, levels = orderSamples)
  ) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = 'Biopsy site', fill = biopsySite.Generalized)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = c('Liver' = '#FF3500', 'Prostate' = '#EDAEAE', 'Bone' = '#FEFEFE', 'Other' = '#4CA947', 'Lymph node' = '#0A6C94'), guide = ggplot2::guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
  themeAnno_Job


## Track - ERG ----

tracks.landscape$ERG <- AbiEnza.Metadata %>%
  dplyr::select(sampleId, hasGenomicERG) %>%
  dplyr::mutate(
    sampleId = factor(sampleId, levels = orderSamples)
  ) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = 'ERG Fusion', fill = hasGenomicERG)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = c('Yes' = '#2F385E', 'No' = 'white')) +
  themeAnno_Job + ggplot2::theme(legend.position = 'none')


## Track - Matching RNA. ----

tracks.landscape$RNA <- AbiEnza.Metadata %>%
  dplyr::select(sampleId, hasMatchingRNA) %>%
  dplyr::mutate(
    sampleId = factor(sampleId, levels = orderSamples)
  ) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = 'Matching RNA', fill = hasMatchingRNA)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = c('Yes' = '#5BA4B8', 'No' = 'white')) +
  themeAnno_Job + ggplot2::theme(legend.position = 'none')

## Combine tracks. ----

tracks.landscape$responderCategory +
  tracks.landscape$treatmentDuration + 
  tracks.landscape$TMB +
  tracks.landscape$SV +
  tracks.landscape$SV.Cat +
  tracks.landscape$Ploidy +
  tracks.landscape$MutSigs.SNV +
  tracks.landscape$MSI +
  tracks.landscape$HRD +
  tracks.landscape$Chromothripsis +
  tracks.landscape$biopsySite +
  tracks.landscape$ERG +
  tracks.landscape$RNA +
  patchwork::plot_layout(guides = 'collect', ncol = 1, heights = c(.25, rep(1, 6), rep(.25, 6))) +
  patchwork::plot_annotation(tag_levels = 'a')


# Differences in genomic landscape ----------------------------------------

## Perform statistical tests. ----
statTests <- list()

statTests$TMB <- AbiEnza.Results$mutationalBurden %>% 
  dplyr::select(sampleId = sample, Genome.TMB) %>% 
  dplyr::inner_join(includedGroups) %>% 
  rstatix::pairwise_wilcox_test(Genome.TMB ~ responderCategory, exact = T, p.adjust.method = 'none', detailed = T, alternative = 'two.sided', paired = F) %>%
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))

statTests$SV <- AbiEnza.Results$mutationalBurden %>% 
  dplyr::select(sampleId = sample, totalSV) %>% 
  dplyr::inner_join(includedGroups) %>% 
  rstatix::pairwise_wilcox_test(totalSV ~ responderCategory, exact = T, p.adjust.method = 'none', detailed = T, alternative = 'two.sided', paired = F) %>%
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))

## Generate figures. ----

plotsDiff <- list()

### Genome-wide TMB per Responder group. ----
plotsDiff$TMB <- AbiEnza.Results$mutationalBurden %>% 
  dplyr::select(sampleId = sample, Genome.TMB) %>% 
  dplyr::inner_join(includedGroups) %>% 
  dplyr::group_by(responderCategory) %>%
  dplyr::mutate(median = round(median(Genome.TMB), 1)) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(., ggplot2::aes(x = responderCategory, y = Genome.TMB, fill = responderCategory, label = median)) +
  gghalves::geom_half_boxplot(side = 'l', alpha = .3, outlier.shape = NA, notch = F, show.legend = F) +
  gghalves::geom_half_point_panel(side = 'r', position = ggbeeswarm::position_quasirandom(width = .15), size = 1.5, shape = 21, color = 'black', lwd = .01) +
  ggplot2::scale_y_continuous(trans = scales::pseudo_log_trans(), expand = c(0,0), breaks = c(0, 2.5, 5, 10, 25, 50, 100), limits = c(0, 175)) +
  ggplot2::scale_fill_manual(values = c('Poor Responder (≤100 days)' = '#E69F00', 'Good Responder (>100 days)' = '#019E73'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
  ggplot2::stat_summary(fun = median, colour = 'black', geom = 'text', size = 3, show.legend = FALSE, vjust = -4, angle = 90) +
  ggplot2::labs(x = NULL, y = 'Tumor Mutational Burden<br><span style = "font-size:5pt">(Genome-wide; log<sub>10</sub>)</span>') +
  ggpubr::geom_bracket(aes(xmin = group1, xmax = group2, label = p.adj.signif, fill = NULL, color = NULL, shape = NULL), data = statTests$TMB, y.position = 4.7, step.increase = .02, tip.length = .01) +
  theme_Job

### Total nr. of SV per Responder group. ----
plotsDiff$SV <- AbiEnza.Results$mutationalBurden %>% 
  dplyr::select(sampleId = sample, totalSV) %>% 
  dplyr::inner_join(includedGroups) %>% 
  dplyr::group_by(responderCategory) %>%
  dplyr::mutate(median = round(median(totalSV), 0)) %>%
  dplyr::ungroup() %>%
  ggplot2::ggplot(., ggplot2::aes(x = responderCategory, y = totalSV, fill = responderCategory, label = median)) +
  gghalves::geom_half_boxplot(side = 'l', alpha = .3, outlier.shape = NA, notch = F, show.legend = F) +
  gghalves::geom_half_point_panel(side = 'r', position = ggbeeswarm::position_quasirandom(width = .15), size = 1.5, shape = 21, color = 'black', lwd = .01) +
  ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(0, 100, 250, 500, 1000, 1500, 2000), limits = c(0, 2250)) +
  ggplot2::scale_fill_manual(values = c('Poor Responder (≤100 days)' = '#E69F00', 'Good Responder (>100 days)' = '#019E73'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, nrow = 1, keywidth = 0.5, keyheight = 0.5)) +
  ggplot2::stat_summary(fun = median, colour = 'black', geom = 'text', size = 3, show.legend = FALSE, vjust = -4, angle = 90) +
  ggplot2::labs(x = NULL, y = 'Nr. of structural variants<br><span style = "font-size:5pt">(Genome-wide)') +
  ggpubr::geom_bracket(aes(xmin = group1, xmax = group2, label = p.adj.signif, fill = NULL, color = NULL, shape = NULL), data = statTests$SV, y.position = 2100, step.increase = .02, tip.length = .01) +
  theme_Job


### Combine plots. ----

plotsDiff$TMB + plotsDiff$SV +   
  patchwork::plot_layout(guides = 'collect', ncol = 2) +
  patchwork::plot_annotation(tag_levels = 'a')
