# Author:    Job van Riet
# Date:      06-01-22
# Function:  Performs analysis on the repeat-biopsies.

# Libraries ----

library(R2CPCT)
library(DESeq2)

# Retrieve RNA-Seq data. ----

# Load metadata of the Abi/Enza-treated patients.
load('/mnt/onco0002/repository/HMF/DR71/Dec2021/RData/AbiEnza.Metadata.RData')

# Load metadata of the DR-071 cohort.
load('/mnt/onco0002/repository/HMF/DR71/Dec2021/RData/DR71.MetaData.RData')

load('/mnt/onco0002/repository/HMF/DR71/Dec2021/RData/AbiEnzaRNA.RData')

# Determine repeated biopsies.
repeatedBiopsies <- DR71.MetaData$sampleInfo %>% 
  dplyr::filter(
    hmfPatientId %in% AbiEnza.Metadata$hmfPatientId, 
    !hmfSampleId %in% AbiEnza.Metadata$hmfSampleId, 
    hasMatchingRNA == 'Yes'
  )

# Previous biopsies.
# c('DRUP01050020T', 'DRUP01070071T', 'DRUP01070057T', 'DRUP01010097T', 'CPCT02070055TII', 'CPCT02020351TII', 'CPCT02070107TII', 'CPCT02010692TII', 'CPCT02140041TII', 'CPCT02140041TIII')
