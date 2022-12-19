# Author:    Job van Riet
# Date:      27-07-22
# Function:  Overview of the predictions.

# Libraries ----

library(R2CPCT)
library(ggplot2)
library(showtext)
library(patchwork)

# Load themes.
source('R/3.Figures/misc_Themes.R')

# Load metadata of the Abi/Enza-treated patients.
load('~/Downloads/AbiEnza.Metadata.RData')

# Retrieve WGS-data.
load('~/Downloads/AbiEnza.Results.RData')

# Read prediction metrics.
dataPred <- readxl::read_xlsx('Misc/Suppl. Table 1 - OverviewOfData.xls', sheet = 'Predictions', skip = 1, trim_ws = T) %>% 
    dplyr::filter(subgroupCohort != 'Training') %>% 
    dplyr::inner_join(AbiEnza.Metadata, by = 'hmfSampleId') %>% 
    dplyr::mutate(Prediction_Clinicogenomics_Ambi = dplyr::if_else(Clinicogenomics_probability_Bad >= .6 | Clinicogenomics_probability_Good >= .6, Prediction_Clinicogenomics, 'Ambiguous Prediction'))


# Figure ----

tracks.landscape <- list()

# Order samples on predictive value (Bad response)
orderSamples <- dataPred %>% dplyr::arrange(-Clinicogenomics_probability_Bad) %>% dplyr::pull(sampleId)

## Track - Responder group ----

tracks.landscape$Responder <- dataPred %>%
  dplyr::mutate(sampleId = factor(sampleId, levels = orderSamples)) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = 'Responder class', fill = Responder)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(dataPred$Responder), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  themeAnno_Job

tracks.landscape$Responder.Predicted <- dataPred %>%
  dplyr::mutate(sampleId = factor(sampleId, levels = orderSamples)) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = 'Predicted class', fill = Prediction_Clinicogenomics_Ambi)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(dataPred$Prediction_Clinicogenomics_Ambi), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  themeAnno_Job


## Prior therapies ----
tracks.landscape$priorAbiEnza <- dataPred %>%
  dplyr::mutate(sampleId = factor(sampleId, levels = orderSamples)) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = 'Prior Abi./Enza.', fill = priorAbiEnza)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T, ) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = c('Yes' = 'black'), na.value = 'white', guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  themeAnno_Job

tracks.landscape$priorChemo <- dataPred %>%
  dplyr::mutate(sampleId = factor(sampleId, levels = orderSamples)) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = 'Prior Doc./Caba.', fill = priorChemo)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T, ) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = c('Yes' = 'black'), na.value = 'white', guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  themeAnno_Job

## Track - Predictive value ----

tracks.landscape$Prediction <- dataPred %>%
  dplyr::mutate(sampleId = factor(sampleId, levels = orderSamples)) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = Clinicogenomics_probability_Good - .5, fill = Clinicogenomics_probability_Good - .5)) +
  ggplot2::geom_bar(stat = 'identity', lwd = .1, color = 'black', width = .8) +
  ggplot2::scale_y_continuous(limits = c(-.5, .5), expand = c(0,0), labels = c('100% - Poor responder', '75% - Poor responder', 'No predictive distinction', '75% - Good responder', '100% - Good responder')) +
  ggplot2::scale_fill_gradient2('Prediction', low = "#8B0000", mid = 'white', high = "#008080", midpoint = 0, limits = c(-.5, .5), labels = c('100% - Poor responder', '75% - Poor responder', 'No predictive distinction', '75% - Good responder', '100% - Good responder'), guide = ggplot2::guide_colorbar(direction = 'vertical', title.position = 'top', title.hjust = 0.5, barwidth = 1, barheight = 5)) +
  ggplot2::labs(y = 'Predictive probability', x = NULL) +
  themeTrack_Job + theme(axis.text.x = ggplot2::element_blank())


## Track - Treatment duration ----

tracks.landscape$treatmentDuration <- dataPred %>%
  dplyr::distinct(sampleId, Treatment, treatmentduration_days) %>%
  dplyr::mutate(sampleId = factor(sampleId, levels = orderSamples)) %>%
  #Plot.
  ggplot2::ggplot(., ggplot2::aes(x = sampleId, xend = sampleId, yend = 0, y = treatmentduration_days, fill = Treatment)) +
  ggplot2::geom_segment() +
  ggplot2::geom_point(shape = 21, size = 1.25) +
  ggplot2::scale_y_continuous(expand = c(0,0), trans = scales::sqrt_trans(), breaks = c(0, 100, 250, 500, 1000, 2000, 2500), limits = c(0, 2501)) +
  ggplot2::geom_hline(yintercept = 180, color = 'black', lty = 15, lwd = 1) +
  ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(dataPred$Treatment), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  ggplot2::labs(y = 'Treatment duration<br><span style = "font-size:5pt">(in days; √)</span>') +
  themeTrack_Job


## Track - Genome-wide TMB ----

tracks.landscape$TMB <- AbiEnza.Results$mutationalBurden %>%
  dplyr::distinct(sample, Genome.TMB) %>%
  dplyr::filter(sample %in% orderSamples) %>%
  dplyr::mutate(sample = factor(sample, levels = orderSamples)) %>%
  #Plot.
  ggplot2::ggplot(., ggplot2::aes(x = sample, xend = sample, yend = 0, y = Genome.TMB)) +
  ggplot2::geom_segment() +
  ggplot2::geom_point(shape = 21, fill = 'red', size = 1.25) +
  ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(0, 5, 10, 25, 50, 75), limits = c(0, 75)) +
  ggplot2::geom_hline(yintercept = 10, color = 'black', lty = 'dotted', lwd = .5) +
  ggplot2::labs(y = 'Tumor Mutational Burden<br><span style = "font-size:5pt">(Genome-wide; √)</span>') +
  ggplot2::coord_trans(y = scales::sqrt_trans()) +
  themeTrack_Job

## Track - Overview of structural variants. ----

tracks.landscape$SV <- AbiEnza.Results$mutationalBurden %>%
  dplyr::distinct(sample, totalSV) %>%
  dplyr::filter(sample %in% orderSamples) %>%
  dplyr::mutate(sample = factor(sample, levels = orderSamples)) %>%
  #Plot.
  ggplot2::ggplot(., ggplot2::aes(x = sample, xend = sample, yend = 0, y = totalSV)) +
  ggplot2::geom_segment() +
  ggplot2::geom_point(shape = 21, fill = '#00AF66', size = 1.25) +
  ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100, 250, 500, 1000, 2000), limits = c(0, 2000)) +
  ggplot2::labs(y = 'No. of structural variants<br><span style = "font-size:5pt">(Genome-wide; √)') +
  ggplot2::coord_trans(y = scales::sqrt_trans()) +
  themeTrack_Job

tracks.landscape$SV.DEL <- AbiEnza.Results$mutationalBurden %>%
  dplyr::distinct(sample, SV.DEL) %>%
  dplyr::filter(sample %in% orderSamples) %>%
  dplyr::mutate(sample = factor(sample, levels = orderSamples)) %>%
  #Plot.
  ggplot2::ggplot(., ggplot2::aes(x = sample, xend = sample, yend = 0, y = SV.DEL)) +
  ggplot2::geom_segment() +
  ggplot2::geom_point(shape = 21, fill = '#ff8c00', size = 1.25) +
  ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100, 250, 500, 1000), limits = c(0, 1000)) +
  ggplot2::labs(y = 'No. of Deletions<br><span style = "font-size:5pt">(Genome-wide; √)') +
  ggplot2::coord_trans(y = scales::sqrt_trans()) +
  themeTrack_Job

tracks.landscape$SV.DUP <- AbiEnza.Results$mutationalBurden %>%
  dplyr::distinct(sample, SV.DUP) %>%
  dplyr::filter(sample %in% orderSamples) %>%
  dplyr::mutate(sample = factor(sample, levels = orderSamples)) %>%
  #Plot.
  ggplot2::ggplot(., ggplot2::aes(x = sample, xend = sample, yend = 0, y = SV.DUP)) +
  ggplot2::geom_segment() +
  ggplot2::geom_point(shape = 21, fill = '#fc6769', size = 1.25) +
  ggplot2::scale_y_continuous(expand = c(0,0), breaks = c(0, 50, 100, 250, 500, 1000), limits = c(0, 1000)) +
  ggplot2::labs(y = 'No. of Tandem duplications<br><span style = "font-size:5pt">(Genome-wide; √)') +
  ggplot2::coord_trans(y = scales::sqrt_trans()) +
  themeTrack_Job


## Track - MSI ----

tracks.landscape$MSI <- AbiEnza.Results$mutationalBurden %>%
  dplyr::select(sample, msStatus) %>%
  dplyr::filter(sample %in% orderSamples) %>%
  dplyr::mutate(
    sample = factor(sample, levels = orderSamples)
  ) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sample, y = 'MSI', fill = msStatus)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = c('MSI' = '#f9b320', 'MSS' = 'white')) +
  themeAnno_Job + ggplot2::theme(legend.position = 'none')


## Track - CHORD. ----

tracks.landscape$HRD <- AbiEnza.Results$mutationalBurden %>%
  dplyr::select(sample, hr_status) %>%
  dplyr::filter(sample %in% orderSamples) %>%
  dplyr::mutate(
    sample = factor(sample, levels = orderSamples),
    isHRD = ifelse(hr_status == 'HR_deficient', 'HR-deficient', 'HR-proficient')
  ) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sample, y = 'HRD', fill = isHRD)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = c('HR-deficient' = '#0079c2', 'HR-proficient' = 'white')) +
  themeAnno_Job + ggplot2::theme(legend.position = 'none')


## Track - Chromothripsis. ----

tracks.landscape$Chromothripsis <- AbiEnza.Results$mutationalBurden %>%
  dplyr::filter(sample %in% orderSamples) %>%
  dplyr::select(sample, hasChromothripsis) %>%
  dplyr::mutate(
    sample = factor(sample, levels = orderSamples)
  ) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sample, y = 'Chromothripsis', fill = hasChromothripsis)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = c('Yes' = '#835244', 'No' = 'white')) +
  themeAnno_Job + ggplot2::theme(legend.position = 'none')


## Track - Biopsy site. ----

tracks.landscape$biopsySite <- dataPred %>%
  dplyr::select(sampleId, Biopsysite_consolidated) %>%
  dplyr::mutate(
    sampleId = factor(sampleId, levels = orderSamples)
  ) %>%
  ggplot2::ggplot(., ggplot2::aes(x = sampleId, y = 'Biopsy site', fill = Biopsysite_consolidated)) +
  ggplot2::geom_tile(width = .8, colour = 'grey25', lwd = .25, na.rm = T) +
  ggplot2::labs(y = NULL, x = NULL) +
  ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(dataPred$Biopsysite_consolidated), guide = guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 2, keywidth = 0.5, keyheight = 0.5)) +
  themeAnno_Job


## Combine tracks. ----

svglite::svglite(file = 'Fig6.svg', width = 11, height = 8.5)
tracks.landscape$Prediction + 
  tracks.landscape$Responder.Predicted +
  tracks.landscape$Responder +
  tracks.landscape$treatmentDuration + 
  tracks.landscape$TMB +
  tracks.landscape$SV +
  tracks.landscape$SV.DEL +
  tracks.landscape$SV.DUP +
  tracks.landscape$priorAbiEnza +
  tracks.landscape$priorChemo +
  tracks.landscape$MSI +
  tracks.landscape$HRD +
  tracks.landscape$Chromothripsis +
  tracks.landscape$biopsySite +
  patchwork::plot_layout(guides = 'collect', ncol = 1, heights = c(1.2, .15, .15, 1, rep(1, 4), rep(.15, 6))) +
  patchwork::plot_annotation(tag_levels = 'a')
dev.off()
