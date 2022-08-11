# Author:    Job van Riet
# Function:
# Generate Kaplan-Meijers based on the prediction classes.

# Libraries and data. ----

library(dplyr)
library(ggplot2)
library(survminer)
library(gtsummary)
library(extrafont)
library(patchwork)

# Load ggplot2 themes.
source('R/3.Figures/misc_Themes.R')
source('R/3.Figures/misc_Functions.R')


# Load metadata of the Abi/Enza-treated patients.
load('~/Downloads/AbiEnza.Metadata.RData')


# Read prediction metrics.
dataPred <- readxl::read_xlsx('Misc/Suppl. Table 1 - OverviewOfData.xls', sheet = 'Predictions', skip = 1, trim_ws = T) %>% 
    dplyr::filter(subgroupCohort != 'Training') %>% 
    dplyr::inner_join(AbiEnza.Metadata, by = 'hmfSampleId') %>% 
    dplyr::mutate(
        Prediction_Transcriptomics = dplyr::if_else(Transcriptomics_probability_Bad >= .6, 'Predicted - Bad Responder', dplyr::if_else(Transcriptomics_probability_Good >= .6, 'Predicted - Good Responder', 'Ambiguous Prediction')),
        Prediction_Genomics = dplyr::if_else(WGS_probability_Bad >= .6, 'Predicted - Bad Responder', dplyr::if_else(WGS_probability_Good >= .6, 'Predicted - Good Responder', 'Ambiguous Prediction')),
        Prediction_Clinicogenomics_Ambi = dplyr::if_else(Clinicogenomics_probability_Bad >= .6 | Clinicogenomics_probability_Good >= .6, Prediction_Clinicogenomics, 'Ambiguous Prediction'),
        treatmentStatus = dplyr::if_else(Treatment_Ongoing == 'Yes', 0, 1)
    )


# Survival Analysis (Cox regression) ----

plotFits <- list()

## ClinicoGenomics (Treatment duration) ----
plotFits$clinicogenomics.Two <- plotKM.Treatment(
    survminer::surv_fit(formula = survival::Surv(treatmentduration_days, treatmentStatus) ~ Prediction_Clinicogenomics, data = dataPred),
    ylim = 2550, palette = c( '#8B0000', '#008080')
)

plotFits$clinicogenomics.Three <- plotKM.Treatment(
    survminer::surv_fit(formula = survival::Surv(treatmentduration_days, treatmentStatus) ~ Prediction_Clinicogenomics_Ambi, data = dataPred),
    ylim = 2550, palette = c('orange', '#8B0000', '#008080')
)

svglite::svglite(file = 'Fig6_Kaplan.svg', width = 5, height = 8.5)
plotFits$clinicogenomics.Two$plot + plotFits$clinicogenomics.Two$table +
    plotFits$clinicogenomics.Three$plot + plotFits$clinicogenomics.Three$table + 
    patchwork::plot_layout(ncol = 1, heights = c(1, .15, 1, .25))
dev.off()


## ClinicoGenomics (Subgroups) ----
plotFits$clinicogenomics.Group1 <- plotKM.Treatment(
    survminer::surv_fit(formula = survival::Surv(treatmentduration_days, treatmentStatus) ~ Prediction_Clinicogenomics_Ambi, data = dataPred %>% dplyr::filter(Number_prior_treatment_lines.x <= 1)),
    ylim = 2550, palette = c('orange', '#8B0000', '#008080')
)

plotFits$clinicogenomics.Group2 <- plotKM.Treatment(
    survminer::surv_fit(formula = survival::Surv(treatmentduration_days, treatmentStatus) ~ Prediction_Clinicogenomics_Ambi, data = dataPred %>% dplyr::filter(Number_prior_treatment_lines.x >= 2)),
    ylim = 2550, palette = c('orange', '#8B0000', '#008080')
)

plotFits$clinicogenomics.Group3 <- plotKM.Treatment(
    survminer::surv_fit(formula = survival::Surv(treatmentduration_days, treatmentStatus) ~ Prediction_Clinicogenomics_Ambi, data = dataPred %>% dplyr::filter(Prior_Enzalutamide.x == 'No')),
    ylim = 2550, palette = c('orange', '#8B0000', '#008080')
)

plotFits$clinicogenomics.Group4 <- plotKM.Treatment(
    survminer::surv_fit(formula = survival::Surv(treatmentduration_days, treatmentStatus) ~ Prediction_Clinicogenomics_Ambi, data = dataPred %>% dplyr::filter(Prior_Enzalutamide.x == 'Yes')),
    ylim = 2550, palette = c('orange', '#8B0000', '#008080')
)

layout = 'AB
EF
CD
GH'

svglite::svglite(file = 'SupplFig5.svg', width = 10, height = 8.5)
plotFits$clinicogenomics.Group1$plot + plotFits$clinicogenomics.Group2$plot + plotFits$clinicogenomics.Group3$plot + plotFits$clinicogenomics.Group4$plot +
    plotFits$clinicogenomics.Group1$table + plotFits$clinicogenomics.Group2$table + plotFits$clinicogenomics.Group3$table + plotFits$clinicogenomics.Group4$table +
    patchwork::plot_layout(design = layout, heights = c(1, .2, 1, .2))
dev.off()


## WGS/WTS-only models. ----

plotFits$WGS.Internal <- plotKM.Treatment(
    survminer::surv_fit(formula = survival::Surv(treatmentduration_days, treatmentStatus) ~ Prediction_Genomics, data = dataPred),
    ylim = 2550, palette = c('orange', '#8B0000', '#008080')
)

plotFits$WTS.Internal <- plotKM.Treatment(
    survminer::surv_fit(formula = survival::Surv(treatmentduration_days, treatmentStatus) ~ Prediction_Transcriptomics, data = dataPred %>% dplyr::filter(!is.na(dataPred$Transcriptomics_probability_Bad))),
    ylim = 2550, palette = c('orange', '#8B0000', '#008080')
)

svglite::svglite(file = 'SupplFig6.svg', width = 9.5, height = 4.5)
plotFits$WGS.Internal$plot + plotFits$WTS.Internal$plot +
    plotFits$WGS.Internal$table + plotFits$WTS.Internal$table + 
    patchwork::plot_layout(ncol = 2, heights = c(1, .25))
dev.off()


## ClinicoGenomics (OS) ----

dataOS <- readxl::read_xlsx('Misc/Suppl. Table 1 - OverviewOfData.xls', sheet = 'Predictions', skip = 1, trim_ws = T) %>% 
    dplyr::inner_join(AbiEnza.Metadata, by = 'hmfSampleId') %>% 
    dplyr::mutate(
        Prediction_Clinicogenomics_Ambi = dplyr::if_else(Clinicogenomics_probability_Bad >= .6 | Clinicogenomics_probability_Good >= .6, Prediction_Clinicogenomics, 'Ambiguous Prediction')
    )

plotFits$clinicogenomics.OS.All.True <- plotKM.OS(
    survminer::surv_fit(formula = survival::Surv(survTimeInDays, Death) ~ Responder, data = dataOS),
    ylim = 2550, palette = c('orange', '#8B0000', '#008080')
)

plotFits$clinicogenomics.OS.Training.Pred <- plotKM.OS(
    survminer::surv_fit(formula = survival::Surv(survTimeInDays, Death) ~ Prediction_Clinicogenomics_Ambi, data = dataOS %>% dplyr::filter(subgroupCohort == 'Training')),
    ylim = 2550, palette = c('orange', '#8B0000', '#008080')
)

plotFits$clinicogenomics.OS.Validation.Pred <- plotKM.OS(
    survminer::surv_fit(formula = survival::Surv(survTimeInDays, Death) ~ Prediction_Clinicogenomics_Ambi, data = dataOS %>% dplyr::filter(subgroupCohort != 'Training')),
    ylim = 2550, palette = c('orange', '#8B0000', '#008080')
)

svglite::svglite(file = 'SupplFig7.svg', width = 16, height = 5)
layout = 'ABC
DEF'
plotFits$clinicogenomics.OS.All.True$plot + plotFits$clinicogenomics.OS.Training.Pred$plot + plotFits$clinicogenomics.OS.Validation.Pred$plot +
    plotFits$clinicogenomics.OS.All.True$table + plotFits$clinicogenomics.OS.Training.Pred$table + plotFits$clinicogenomics.OS.Validation.Pred$table +
    patchwork::plot_layout(design = layout, heights = c(1, .2)) + patchwork::plot_annotation(tag_levels = 'a')
dev.off()
