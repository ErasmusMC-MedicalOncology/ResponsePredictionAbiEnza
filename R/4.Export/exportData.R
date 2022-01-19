# Author:    Job van Riet
# Date:      19-01-22
# Function:  Export results to Excel.

# Libraries ----

library(R2CPCT)

# Import data.
load('/mnt/onco0002/repository/HMF/DR71/Dec2021/RData/AbiEnza.Metadata.RData')
load('/mnt/onco0002/repository/HMF/DR71/Dec2021/RData/AbiEnza.Results.RData')
load('/mnt/onco0002/repository/HMF/DR71/Dec2021/RData/AbiEnza.RNASeq.RData')


# Export results --------------------------------------------------------------------------------------------------

# Open workbook.
wb <- openxlsx::createWorkbook(creator = 'Job van Riet', subject = 'AbiEnza', title = 'AbiEnza', category = 'CPCT-02')

## Sample metadata ----
openxlsx::addWorksheet(wb, 'Sample information')
data.Meta <- AbiEnza.Metadata %>% 
  dplyr::mutate(sampleId = NULL, patientIdentifier= NULL, setName = NULL, startDate_treatment = NULL, endDate_treatment = NULL, Treatment_ongoing_19122021 = NULL, Time_between_biopsy_tm_days = NULL, biopsyDate = NULL, biopsySite = NULL, biopsyLocation = NULL) %>% 
  dplyr::mutate(biopsySite = Biopsysite_consolidated, Biopsysite_consolidated = NULL)

openxlsx::writeDataTable(wb, sheet = 'Sample information', x = data.Meta)

## Mutation burden ----
openxlsx::addWorksheet(wb, 'Mutation Burden')
data.TMB <- AbiEnza.Results$mutationalBurden %>% 
  dplyr::inner_join(AbiEnza.Metadata %>% dplyr::select(sample = sampleId, hmfSampleId)) %>%
  dplyr::mutate(sample = hmfSampleId, hmfSampleId = NULL)

openxlsx::writeDataTable(wb, sheet = 'Mutation Burden', x = data.TMB)

## dNdS ----
openxlsx::addWorksheet(wb, 'dNdS')
data.dNdS <- dplyr::bind_rows(
  AbiEnza.Results$dNdS$finalOutput %>% dplyr::filter(qmis_loc <= .1 | qmis_loc<= .1 | qall_loc <= .1  | qmis_cv <= .1 | qallsubs_cv <= .1 | qglobal_cv <= .1) %>% dplyr::mutate(Cohort = 'Entire cohort'),
  AbiEnza.Results$dNdS.GoodResponders$finalOutput %>% dplyr::filter(qmis_loc <= .1 | qmis_loc<= .1 | qall_loc <= .1  | qmis_cv <= .1 | qallsubs_cv <= .1 | qglobal_cv <= .1) %>% dplyr::mutate(Cohort = 'Good Responders'),
  AbiEnza.Results$dNdS.BadResponders$finalOutput %>% dplyr::filter(qmis_loc <= .1 | qmis_loc<= .1 | qall_loc <= .1  | qmis_cv <= .1 | qallsubs_cv <= .1 | qglobal_cv <= .1) %>% dplyr::mutate(Cohort = 'Bad Responders')
) %>% 
  dplyr::select(SYMBOL, ENSEMBL, Cohort, dplyr::everything())

openxlsx::writeDataTable(wb, sheet = 'dNdS', x = data.dNdS)

## GISTIC2 ----
openxlsx::addWorksheet(wb, 'GISTIC2')

data.GISTIC <- dplyr::bind_rows(
  data.frame(AbiEnza.Results$GISTIC2.AllSamples$gisticNarrowPeaksWithAnno) %>% dplyr::select(Unique.Name, Descriptor, Peak.Limits, Region.Limits, q.values, Residual.q.values.after.removing.segments.shared.with.higher.peaks, nGenes.GENCODE, overlapGenes.Drivers) %>% dplyr::mutate(Cohort = 'All samples'),
  data.frame(AbiEnza.Results$GISTIC2.GoodResponders$gisticNarrowPeaksWithAnno) %>% dplyr::select(Unique.Name, Descriptor, Peak.Limits, Region.Limits, q.values, Residual.q.values.after.removing.segments.shared.with.higher.peaks, nGenes.GENCODE, overlapGenes.Drivers) %>% dplyr::mutate(Cohort = 'Good responders'),
  data.frame(AbiEnza.Results$GISTIC2.BadResponders$gisticNarrowPeaksWithAnno) %>% dplyr::select(Unique.Name, Descriptor, Peak.Limits, Region.Limits, q.values, Residual.q.values.after.removing.segments.shared.with.higher.peaks, nGenes.GENCODE, overlapGenes.Drivers) %>% dplyr::mutate(Cohort = 'Bad responders')
)

openxlsx::writeDataTable(wb, sheet = 'GISTIC2', x = data.GISTIC)

## Mutational Signatures ----
openxlsx::addWorksheet(wb, 'Mutational signatures')

data.MutSigs <- AbiEnza.Results$mutSigs$SNV$relativeContribution %>% reshape2::dcast(sampleId ~ Signature, value.var = 'relContribution') %>%
  dplyr::inner_join(AbiEnza.Results$mutSigs$InDel$relativeContribution %>% reshape2::dcast(sampleId ~ Signature, value.var = 'relContribution')) %>%
  dplyr::inner_join(AbiEnza.Results$mutSigs$DBS$relativeContribution %>% reshape2::dcast(sampleId ~ Signature, value.var = 'relContribution')) %>%
  dplyr::inner_join(AbiEnza.Metadata %>% dplyr::select(sampleId, hmfSampleId)) %>%
  dplyr::mutate(sampleId = hmfSampleId, hmfSampleId = NULL) %>% 
  dplyr::select(sampleId, dplyr::everything())

openxlsx::writeDataTable(wb, sheet = 'Mutational signatures', x = data.MutSigs)

## Combined Report ----
openxlsx::addWorksheet(wb, 'Mutation Report')
data.Report <- AbiEnza.Results$combinedReport %>% 
  dplyr::filter(isMutant) %>% 
  dplyr::inner_join(AbiEnza.Metadata %>% dplyr::select(sample = sampleId, hmfSampleId, Responder)) %>%
  dplyr::mutate(sample = hmfSampleId, hmfSampleId = NULL, Peak.Amplitude = NULL, isMutant = NULL, chr = NULL, geneId = NULL) %>% 
  dplyr::select(sample, Responder, SYMBOL, ENSEMBL, dplyr::everything())

openxlsx::writeDataTable(wb, sheet = 'Mutation Report', x = data.Report)

## DESeq2 ----

openxlsx::addWorksheet(wb, sheet = 'Results - DESeq2')
data.DESeq2 <- AbiEnza.RNASeq$DESeq2Results %>% 
  dplyr::filter(isSig == 'Significant') %>% 
  dplyr::mutate(contrast = 'Bad vs. Good responders', isSig = NULL) %>% 
  dplyr::select(SYMBOL, ENSEMBL, dplyr::everything())

openxlsx::writeDataTable(wb, sheet = 'Results - DESeq2', x = data.DESeq2)

# Write to file ---------------------------------------------------------------------------------------------------

openxlsx::saveWorkbook(wb, file = '~/test/AbiEnza.xls', overwrite = T)