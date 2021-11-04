# Author:    Job van Riet
# Date:      4-11-21
# Function:  Export results to Excel.

# Libraries ----

library(R2CPCT)

# Import data.
load('/mnt/onco0002/repository/HMF/DR71/Oct2021/RData/AbiEnza.Metadata.RData')
load('/mnt/onco0002/repository/HMF/DR71/Oct2021/RData/AbiEnza.Results.RData')


# Export results --------------------------------------------------------------------------------------------------

# Open workbook.
wb <- openxlsx::createWorkbook(creator = 'Job van Riet', subject = 'AbiEnza', title = 'AbiEnza', category = 'CPCT-02')

## Sample metadata ----
openxlsx::addWorksheet(wb, 'Sample information')
openxlsx::writeDataTable(wb, sheet = 'Sample information', x = AbiEnza.Metadata)

## Mutation burden ----
openxlsx::addWorksheet(wb, 'Mutation Burden')
data.TMB <- AbiEnza.Results$mutationalBurden %>% 
  dplyr::inner_join(AbiEnza.Metadata %>% dplyr::select(sample = sampleId, hmfSampleId)) %>%
  dplyr::mutate(sample = hmfSampleId, hmfSampleId = NULL)

openxlsx::writeDataTable(wb, sheet = 'Mutation Burden', x = data.TMB)

## dNdS ----
openxlsx::addWorksheet(wb, 'dNdS')
openxlsx::writeDataTable(wb, sheet = 'dNdS', x = AbiEnza.Results$dNdS$finalOutput %>% dplyr::filter(qglobal_cv <= 0.1 | qallsubs_cv <= 0.1))

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
  dplyr::mutate(sampleId = hmfSampleId, hmfSampleId = NULL)

openxlsx::writeDataTable(wb, sheet = 'Mutational signatures', x = data.MutSigs)

## Mut. Excl. ----
openxlsx::addWorksheet(wb, 'Mut. Excl.')
openxlsx::writeDataTable(wb, sheet = 'Mut. Excl.', x = AbiEnza.Results$differencesWGS$mutExcl)

## Combined Report ----
openxlsx::addWorksheet(wb, 'Mutation Report')
data.Report <- AbiEnza.Results$combinedReport %>% 
  dplyr::inner_join(AbiEnza.Metadata %>% dplyr::select(sample = sampleId, hmfSampleId)) %>%
  dplyr::mutate(sample = hmfSampleId, hmfSampleId = NULL)

openxlsx::writeDataTable(wb, sheet = 'Mutation Report', x = data.Report)

# Write to file ---------------------------------------------------------------------------------------------------

openxlsx::saveWorkbook(wb, file = '~/test/AbiEnza.xls', overwrite = T)