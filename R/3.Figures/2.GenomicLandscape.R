# Author:    Job van Riet
# Date:      2-11-21
# Function:  Genomic overview of the Abi/Enza-treated patients.

# Libraries ----

library(R2CPCT)
library(ggplot2)
library(extrafont)
library(patchwork)

# Load themes.
source('R/3.Figures/misc_Themes.R')

# Load metadata of the Abi/Enza-treated patients.
load('/mnt/onco0002/repository/HMF/DR71/Oct2021/RData/AbiEnza.Metadata.RData')

# Retrieve WGS-data.
load('/mnt/onco0002/repository/HMF/DR71/Oct2021/RData/AbiEnza.Results.RData')


# Main Figure - Genomic Landscape -----------------------------------------

tracks.landscape <- list()

# Sort on Abi/Enza treatment duration.
orderSamples <- AbiEnza.Metadata %>%
  dplyr::inner_join(AbiEnza.Results$mutationalBurden, by = c('sampleId' = 'sample')) %>% 
  dplyr::mutate(Responder = factor(Responder, levels = c('Good Responder (>100 days)', 'Bad Responder (≤100 days)'))) %>% 
  dplyr::arrange(Responder, -treatmentDurationInDays) %>%
  dplyr::pull(sampleId)


## Track - Responder group ----

tracks.landscape$Responder <- AbiEnza.Metadata %>%
  dplyr::distinct(sampleId, Responder) %>%
  dplyr::mutate(sampleId = factor(sampleId, levels = orderSamples)) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = 'Responder class', fill = Responder)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(AbiEnza.Metadata$Responder), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
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
  ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(AbiEnza.Metadata$treatment), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  ggplot2::labs(y = 'Treatment duration<br><span style = "font-size:5pt">(in days)</span>') +
  themeTrack_Job


## Track - Genome-wide TMB ----

tracks.landscape$TMB <- AbiEnza.Results$mutationalBurden %>%
  dplyr::distinct(sample, Genome.TMB) %>%
  dplyr::mutate(sample = factor(sample, levels = orderSamples)) %>%
  #Plot.
  ggplot2::ggplot(., ggplot2::aes(x = sample, y = Genome.TMB)) +
  ggplot2::geom_bar(stat = 'identity', fill = '#F75658', color = 'black', lwd = .33, width = .8) +
  ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(0, 2.5, 5, 10, 25, 50, 100), labels = c(0, 2.5, 5, 10, 25, 50, 100), limits = c(0, 115)) +
  ggplot2::geom_hline(yintercept = 10, color = 'black', lty = 'dotted', lwd = .5) +
  ggplot2::labs(y = 'Tumor Mutational Burden<br><span style = "font-size:5pt">(Genome-wide; log<sub>10</sub>)</span>') +
  ggplot2::coord_trans(y = scales::pseudo_log_trans()) +
  themeTrack_Job

## Track - Overview of structural variants. ----

tracks.landscape$SV <- AbiEnza.Results$mutationalBurden %>%
  dplyr::distinct(sample, totalSV) %>%
  dplyr::mutate(sample = factor(sample, levels = orderSamples)) %>%
  #Plot.
  ggplot2::ggplot(., ggplot2::aes(x = sample, y = totalSV)) +
  ggplot2::geom_bar(stat = 'identity', fill = '#00AF66', color = 'black', lwd = .33, width = .8) +
  ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(0, 100, 250, 500, 1000, 1500, 2000), limits = c(0, 2100)) +
  ggplot2::labs(y = 'No. of structural variants<br><span style = "font-size:5pt">(Genome-wide)') +
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
  ggplot2::scale_fill_manual(values = colorPalette, breaks = c('Deletions', 'Inversions', 'Tandem Duplications', 'Translocations', 'Insertions', 'Singles'), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
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

tracks.landscape$MutSigs.SNV <- plotMutationalSignaturesCOSMIC(mutSigs = AbiEnza.Results$mutSigs$SNV, orderSamples = as.character(orderSamples), minContrib = 5, combineSigs = T) + themeTrack_Job


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
  ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(AbiEnza.Metadata$biopsySite.Generalized), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 2, keywidth = 0.5, keyheight = 0.5)) +
  themeAnno_Job


## Track - Matching RNA. ----

tracks.landscape$RNA <- AbiEnza.Metadata %>%
  dplyr::select(sampleId, matchedRNA) %>%
  dplyr::mutate(
    sampleId = factor(sampleId, levels = orderSamples)
  ) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = 'Matching RNA', fill = matchedRNA)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = c('Yes' = '#5BA4B8', 'No' = 'white')) +
  themeAnno_Job + ggplot2::theme(legend.position = 'none')


## Combine tracks. ----

tracks.landscape$Responder +
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
  tracks.landscape$RNA +
  patchwork::plot_layout(guides = 'collect', ncol = 1, heights = c(.2, rep(1, 6), rep(.2, 5))) +
  patchwork::plot_annotation(tag_levels = 'a')

