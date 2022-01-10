# Author:    Job van Riet
# Date:      07-01-22
# Function:  Analyze genomic differences between responders classes.

# Libraries ---------------------------------------------------------------

library(R2CPCT)

# Load metadata of the Abi/Enza-treated patients.
load('/mnt/onco0002/repository/HMF/DR71/Dec2021/RData/AbiEnza.Metadata.RData')

# Retrieve WGS-data.
load('/mnt/onco0002/repository/HMF/DR71/Dec2021/RData/AbiEnza.Results.RData')

# Import required data.
data('driverList', package = 'R2CPCT')

# Add list to contain the results of the differences.
AbiEnza.Results$differencesWGS <- list()

# Differences - Chromosomal arms. ----

AbiEnza.Results$differencesWGS$chromosomalArm <- AbiEnza.Results$GISTIC2.AllSamples$gisticBroadScoresPerArm %>% 
  dplyr::inner_join(AbiEnza.Metadata %>% dplyr::select(sampleId, Responder), by = c('variable' = 'sampleId')) %>% 
  dplyr::filter(!grepl('Ambi', Responder)) %>% 
  dplyr::group_by(`Chromosome Arm`) %>% 
  rstatix::pairwise_wilcox_test(value ~ Responder, p.adjust.method = 'none', detailed = T) %>%
  dplyr::ungroup() %>% 
  rstatix::adjust_pvalue(method = 'BH')


# Differences - Mutational Burdens ----

AbiEnza.Results$differencesWGS$MutationalBurden <- AbiEnza.Results$mutationalBurden %>% 
  dplyr::mutate(SV.SINGLE = NULL) %>% 
  dplyr::select(sampleId = sample, Genome.TMB, totalSV, genomePloidy, dplyr::contains(c('SV.'))) %>% 
  reshape2::melt(id.vars = c('sampleId')) %>% 
  dplyr::inner_join(AbiEnza.Metadata) %>% 
  dplyr::filter(!grepl('Ambi', Responder)) %>% 
  dplyr::group_by(variable) %>% 
  rstatix::pairwise_wilcox_test(value ~ Responder, exact = T, p.adjust.method = 'none', detailed = T, alternative = 'two.sided', paired = F) %>%
  dplyr::ungroup() %>% 
  rstatix::adjust_pvalue(method = 'BH') %>% 
  rstatix::add_significance(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c('***', '**', '*', 'ns'))


# Differences - Mutually-exclusive drivers ----

## Define groups. ----
includedGroups <- AbiEnza.Metadata %>% 
  dplyr::distinct(sampleId, Responder) %>% 
  dplyr::filter(!grepl('Ambi', Responder)) %>% 
  dplyr::group_by(Responder) %>% 
  dplyr::mutate(totalInGroup = dplyr::n_distinct(sampleId)) %>% 
  dplyr::ungroup()

## Retrieve all relevant driver(-like) genes. ----

# Min. aberrations (20% in group)
minSamples.Mut <- AbiEnza.Results$combinedReport %>%
  dplyr::filter(isMutant, ENSEMBL %in% driverList$ENSEMBL) %>% 
  dplyr::inner_join(includedGroups, by = c('sample' = 'sampleId')) %>%
  dplyr::group_by(Responder, ENSEMBL, SYMBOL) %>%
  dplyr::summarise(
    totalMuts = dplyr::n_distinct(sample),
    totalMuts.Rel = totalMuts / totalInGroup
  ) %>%
  dplyr::ungroup() %>%
  dplyr::filter(totalMuts.Rel >= .25) %>%
  dplyr::distinct(ENSEMBL) %>% 
  dplyr::pull(ENSEMBL)

# a priori selected genetic aberrations.
selectedGenes <- c('AR', 'TP53', 'PTEN', 'RB1', 'CTNNB1')

# Retrieve / concatenate all relevant genes.
driverGenes <- unique(c(
  AbiEnza.Results$dNdS.GoodResponders$finalOutput %>% dplyr::filter(qmis_loc <= .05 | qmis_loc<= .05 | qall_loc <= .05  | qmis_cv <= .05 | qallsubs_cv <= .05 | qglobal_cv <= .05) %>% dplyr::pull(ENSEMBL),
  AbiEnza.Results$dNdS.BadResponders$finalOutput %>% dplyr::filter(qmis_loc <= .05 | qmis_loc<= .05 | qall_loc <= .05  | qmis_cv <= .05 | qallsubs_cv <= .05 | qglobal_cv <= .05) %>% dplyr::pull(ENSEMBL),
  minSamples.Mut,
  selectedGenes
))

# Count occurrences. ----

# Count mutant-genes.
mutData <- AbiEnza.Results$combinedReport %>%
  dplyr::filter(ENSEMBL %in% driverGenes | SYMBOL %in% driverGenes) %>% 
  dplyr::inner_join(includedGroups, by = c('sample' = 'sampleId')) %>% 
  dplyr::group_by(Responder, SYMBOL) %>%
  dplyr::summarise(
    totalMut = sum(isMutant),
    noMut = totalInGroup - totalMut
  ) %>%
  dplyr::ungroup() %>% 
  dplyr::filter(!is.na(SYMBOL)) %>%
  dplyr::distinct()

# Complete missing data.
mutData <- mutData %>%
  tidyr::complete(SYMBOL, Responder) %>%
  dplyr::rowwise() %>% 
  dplyr::mutate(
    totalMut = ifelse(is.na(totalMut), 0, totalMut),
    noMut = ifelse(is.na(noMut), unique(includedGroups[includedGroups$Responder == unique(Responder),]$totalInGroup), noMut)
  ) %>% 
  dplyr::ungroup()

# Perform Fisher's Exact Test between responder categories.
fisherData <- do.call(rbind, lapply(unique(mutData$SYMBOL), function(gene){
  
  geneData <- mutData %>%
    dplyr::filter(SYMBOL == gene) %>%
    dplyr::summarise(
      Good.withMut = .data$totalMut[Responder == 'Good Responder (≥180 days)'],
      Good.withoutMut = .data$noMut[Responder == 'Good Responder (≥180 days)'],
      
      Bad.withMut = .data$totalMut[Responder != 'Good Responder (≥180 days)'],
      Bad.withoutMut = .data$noMut[Responder != 'Good Responder (≥180 days)'],
    )
  
  test <- data.frame(
    row.names = c('A', 'B'),
    mut = c(geneData$Good.withMut, geneData$Good.withoutMut),
    noMut = c(geneData$Bad.withMut, geneData$Bad.withoutMut)
  )
  
  geneData$SYMBOL <- gene
  geneData$p <- fisher.test(test, hybrid = F, alternative = 'two', simulate.p.value = T)$p.value
  
  return(geneData)
  
}))

# Check direction and effect size.
fisherData <- fisherData %>% dplyr::mutate(
  effectSize.Good = round((Good.withMut / (Good.withMut + Good.withoutMut)) * 100, 1),
  effectSize.Bad = round((Bad.withMut / (Bad.withMut + Bad.withoutMut)) * 100, 1),
  effectSize = sprintf('%s%% vs. %s%%', effectSize.Good, effectSize.Bad)
)

# Correct for multiple-testing.
fisherData$p.adj <- stats::p.adjust(fisherData$p, method = 'BH')

AbiEnza.Results$differencesWGS$mutExcl <- fisherData %>% dplyr::arrange(p.adj)


# Save to object. ----

save(AbiEnza.Results, file = '/mnt/onco0002/repository/HMF/DR71/Dec2021/RData/AbiEnza.Results.RData')
