# Author:    Job van Riet
# Date:      06-04-21
# Function:  Analyze genomic differences between Good and Poor responders.

# Libraries ---------------------------------------------------------------

library(R2CPCT)

# Load metadata of the Abi/Enza-treated patients.
AbiEnza.Metadata <- readxl::read_xlsx('Misc/Suppl. Table 1 - OverviewOfData.xlsx', sheet = 'Sample overview')

# Retrieve WGS-data.
load('/mnt/data2/hartwig/DR71/Oct2020_AbiEnza/RData/AbiEnza.Results.RData')


# Detect mutually-exclusive aberrations -----------------------------------

# Only check between the major groups.
includedGroups <- AbiEnza.Metadata %>% 
  dplyr::filter(responderCategory != 'Unknown Responder') %>% 
  dplyr::distinct(sampleId, responderCategory) %>% 
  dplyr::group_by(responderCategory) %>% 
  dplyr::mutate(totalInGroup = dplyr::n_distinct(sampleId)) %>% 
  dplyr::ungroup()

# Retrieve all relevant driver(-like) genes. ----

## Min. aberrations (20% in group)
minSamples.Mut <- AbiEnza.Results$combinedReport %>%
  dplyr::filter(isMutant, ENSEMBL %in% R2CPCT::driverList$ENSEMBL) %>% 
  dplyr::inner_join(includedGroups, by = c('sample' = 'sampleId')) %>%
  dplyr::group_by(responderCategory, ENSEMBL, SYMBOL) %>%
  dplyr::summarise(
    totalMuts = dplyr::n_distinct(sample),
    totalMuts.Rel = totalMuts / totalInGroup
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(totalMuts.Rel >= .2) %>%
  dplyr::pull(ENSEMBL)

# Retrieve / concatenate all relevant genes.
driverGenes <- unique(c(
  AbiEnza.Results$dNdS$finalOutput %>% dplyr::filter(qallsubs_cv <= .1 | qglobal_cv <= .1) %>% dplyr::pull(ENSEMBL),
  AbiEnza.Results$dNdS.GoodResponders$finalOutput %>% dplyr::filter(qallsubs_cv <= .1 | qglobal_cv <= .1) %>% dplyr::pull(ENSEMBL),
  AbiEnza.Results$dNdS.PoorResponders$finalOutput %>% dplyr::filter(qallsubs_cv <= .1 | qglobal_cv <= .1) %>% dplyr::pull(ENSEMBL),
  minSamples.Mut
))

# Count occurrences. ----

# Count mutant-genes.
mutData <- AbiEnza.Results$combinedReport %>%
  dplyr::filter(ENSEMBL %in% driverGenes) %>% 
  dplyr::inner_join(includedGroups, by = c('sample' = 'sampleId')) %>% 
  dplyr::group_by(responderCategory, SYMBOL) %>%
  dplyr::summarise(
    totalMut = sum(isMutant),
    noMut = totalInGroup - totalMut
  ) %>%
  dplyr::ungroup()

# Count HRD-positives.
mutHRD <- AbiEnza.Results$mutationalBurden %>%
  dplyr::inner_join(includedGroups, by = c('sample' = 'sampleId')) %>% 
  dplyr::group_by(responderCategory) %>%
  dplyr::summarise(
    totalMut = dplyr::n_distinct(sample[hr_status == 'HR_deficient']),
    noMut = dplyr::n_distinct(sample[hr_status != 'HR_deficient']),
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(SYMBOL = 'HRD')

# Count chromothripsis-positives.
mutChromo <- AbiEnza.Results$mutationalBurden %>%
  dplyr::inner_join(includedGroups, by = c('sample' = 'sampleId')) %>% 
  dplyr::group_by(responderCategory) %>%
  dplyr::summarise(
    totalMut = dplyr::n_distinct(sample[hasChromothripsis == 'Yes']),
    noMut = dplyr::n_distinct(sample[hasChromothripsis != 'Yes']),
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(SYMBOL = 'Chromothripsis')

# Count ERG-positives.
mutERG <- AbiEnza.Metadata %>%
  dplyr::filter(sampleId %in% includedGroups$sampleId) %>% 
  dplyr::group_by(responderCategory) %>%
  dplyr::summarise(
    totalMut = dplyr::n_distinct(sampleId[hasGenomicERG == 'Yes']),
    noMut = dplyr::n_distinct(sampleId[hasGenomicERG != 'Yes']),
  ) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(SYMBOL = 'ERG')

# Combine.
mutData <- rbind(mutData, mutHRD, mutChromo, mutERG) %>% 
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::distinct()

# Complete missing data.
mutData <- mutData %>%
  tidyr::complete(SYMBOL, responderCategory) %>%
  dplyr::mutate(
    totalMut = ifelse(is.na(totalMut), 0, totalMut),
    noMut = ifelse(is.na(noMut), unique(includedGroups[includedGroups$responderCategory == unique(responderCategory),]$totalInGroup), noMut)
  )

# Perform Fisher's Exact Test between responder categories.
fisherData <- do.call(rbind, lapply(unique(mutData$SYMBOL), function(gene){
  
  geneData <- mutData %>%
    dplyr::filter(SYMBOL == gene) %>%
    dplyr::summarise(
      Good.withMut = totalMut[responderCategory == 'Good Responder (>100 days)'],
      Good.withoutMut = noMut[responderCategory == 'Good Responder (>100 days)'],
      
      Poor.withMut = sum(totalMut[responderCategory != 'Good Responder (>100 days)']),
      Poor.withoutMut = sum(noMut[responderCategory != 'Good Responder (>100 days)']),
    )
  
  test <- data.frame(
    row.names = c('A', 'B'),
    mut = c(geneData$Good.withMut, geneData$Good.withoutMut),
    noMut = c(geneData$Poor.withMut, geneData$Poor.withoutMut)
  )
  
  geneData$SYMBOL <- gene
  geneData$p <- fisher.test(test, hybrid = F, alternative = 'two', simulate.p.value = T)$p.value
  
  return(geneData)

}))

# Check direction and effect size.
fisherData <- fisherData %>% dplyr::mutate(
  effectSize.Good = round((Good.withMut / (Good.withMut + Good.withoutMut)) * 100, 1),
  effectSize.Poor = round((Poor.withMut / (Poor.withMut + Poor.withoutMut)) * 100, 1),
  effectSize = sprintf('%s%% vs. %s%%', effectSize.Good, effectSize.Poor)
)

# Correct for multiple-testing.
fisherData$p.adj <- stats::p.adjust(fisherData$p, method = 'BH')


# Differences in Mut. Signatures ------------------------------------------

sigData <- AbiEnza.Results$mutSigs$SNV$relativeContribution %>% 
  dplyr::inner_join(includedGroups) %>% 
  dplyr::group_by(Signature) %>% 
  dplyr::summarise(p = stats::pairwise.wilcox.test(relContribution, responderCategory, p.adjust.method = 'none', paired = FALSE)$p.value) %>% 
  dplyr::ungroup() %>% 
  dplyr::mutate(p.adj = stats::p.adjust(p, method = 'BH'))

