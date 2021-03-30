# Author:    Job van Riet
# Date:      30-03-21
# Function:  Import and processing of the Abi/Enza-treated WGS samples (DR-071).

# Libraries ---------------------------------------------------------------

library(R2CPCT)
library(GenomicRanges)
library(ShatterSeek)

# Initialize a file-logger to store INFO / TRACE messages.
R2CPCT::initializeLogger(file = '~/AbiEnza.log')

# Load metadata of the Abi/Enza-treated patients.
AbiEnza.Metadata <- readxl::read_xlsx('Misc/Suppl. Table 1 - OverviewOfData.xlsx', sheet = 'Sample overview')

# Import / convert WGS samples --------------------------------------------

# Import the WGS data of all samples in the cohort.
AbiEnza.CohortWGS <- R2CPCT::importWGSOfCohort(AbiEnza.Metadata$sampleId, '/mnt/data2/hartwig/DR71/Oct2020/dataHMF/combinedData/', nThreads = 20)
save(AbiEnza.CohortWGS, file = '/mnt/data2/hartwig/DR71/Oct2020_AbiEnza/RData/AbiEnza.CohortWGS.RData')


# Perform GISTIC2 analysis ------------------------------------------------

# Generate the commands to perform GISTIC2 (perform this in your Bash terminal).
# Perform this separately per responder group (Good / Poor).

# All samples.
R2CPCT::performGISTIC2(AbiEnza.CohortWGS$copyNumbers, outputFolder = '/mnt/data2/hartwig/DR71/Oct2020_AbiEnza/results/GISTIC2/AllSamples/')

# Good responders.
R2CPCT::performGISTIC2(AbiEnza.CohortWGS$copyNumbers[GenomicRanges::mcols(AbiEnza.CohortWGS$copyNumbers)$sample %in% (AbiEnza.Metadata %>% dplyr::filter(grepl('Good Responder', responderCategory)) %>% dplyr::pull(sampleId)),], outputFolder = '/mnt/data2/hartwig/DR71/Oct2020_AbiEnza/results/GISTIC2/GoodResponders/')

# Poor responders.
R2CPCT::performGISTIC2(AbiEnza.CohortWGS$copyNumbers[GenomicRanges::mcols(AbiEnza.CohortWGS$copyNumbers)$sample %in% (AbiEnza.Metadata %>% dplyr::filter(grepl('Poor Responder', responderCategory)) %>% dplyr::pull(sampleId)),], outputFolder = '/mnt/data2/hartwig/DR71/Oct2020_AbiEnza/results/GISTIC2/PoorResponders/')

# GISTIC2 command used:
# GISTIC2_2.0.23/gistic2 -b /mnt/data2/hartwig/DR71/Oct2020_AbiEnza/results/GISTIC2/PoorResponders/ -seg /tmp/Rtmp7hzfek/50de7e5bd09e.txt -refgene GISTIC2_2.0.23/refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat -genegistic 1 -gcm extreme -maxseg 4000 -broad 1 -brlen 0.98 -conf 0.95 -rx 0 -cap 3 -saveseg 0 -armpeel 1 -smallmem 0 -res 0.01 -ta 0.3 -td 0.3 -savedata 0 -savegene 1 -qvt 0.1 -twoside 0


# Cohort-wide analysis ----------------------------------------------------

# Initialize a list to contain the results of the WGS analysis.
AbiEnza.Results <- list()

# Generate overview of mutational burden per sample and cohort.
AbiEnza.Results$mutationalBurden <- R2CPCT::determineMutationalBurden(AbiEnza.CohortWGS, minTAF.Muts = 0, minTAF.SV = 0)

# Perform dN/dS on entire cohort and per responder groups.
AbiEnza.Results$dNdS <- R2CPCT::rundNdS(AbiEnza.CohortWGS$somaticVariants)
AbiEnza.Results$dNdS.GoodResponders <- R2CPCT::rundNdS(AbiEnza.CohortWGS$somaticVariants[AbiEnza.CohortWGS$somaticVariants$sample %in% (AbiEnza.Metadata %>% dplyr::filter(grepl('Good Responder', responderCategory)) %>% dplyr::pull(sampleId)),])
AbiEnza.Results$dNdS.PoorResponders <- R2CPCT::rundNdS(AbiEnza.CohortWGS$somaticVariants[AbiEnza.CohortWGS$somaticVariants$sample %in% (AbiEnza.Metadata %>% dplyr::filter(grepl('Poor Responder', responderCategory)) %>% dplyr::pull(sampleId)),])

# Import the GISTIC2-determined recurrent CNA peaks.
AbiEnza.Results$GISTIC2.AllSamples <- R2CPCT::importGISTIC2('/mnt/data2/hartwig/DR71/Oct2020_AbiEnza/results/GISTIC2/AllSamples/')
AbiEnza.Results$GISTIC2.GoodResponders <- R2CPCT::importGISTIC2('/mnt/data2/hartwig/DR71/Oct2020_AbiEnza/results/GISTIC2/GoodResponders/')
AbiEnza.Results$GISTIC2.PoorResponders <- R2CPCT::importGISTIC2('/mnt/data2/hartwig/DR71/Oct2020_AbiEnza/results/GISTIC2/PoorResponders/')

# Generate the gene-level mutational overview.
AbiEnza.Results$combinedReport <- R2CPCT::generateCombinedReport(AbiEnza.CohortWGS, dNdS = AbiEnza.Results$dNdS, GISTIC2 = AbiEnza.Results$GISTIC2.AllSamples, nThreads = 20, mutantsOnly = T)

# Determine mutational signatures.
AbiEnza.Results$mutSigs <- R2CPCT::fitMutSigs(AbiEnza.CohortWGS$somaticVariants, restrictiveFit = F)

# Determine # of Ti/Tv and ratio per sample.
AbiEnza.Results$TiTv <- R2CPCT::determineTiTv(AbiEnza.Results$mutSigs$SNV$mutMatrix)

# Save results.
save(AbiEnza.Results, file = '/mnt/data2/hartwig/DR71/Oct2020_AbiEnza/RData/AbiEnza.Results.RData')


# Close logger ------------------------------------------------------------

ParallelLogger::unregisterLogger()
