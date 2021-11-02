# Author:    Job van Riet
# Date:      27-10-21
# Function:  Import and processing of the Abi/Enza-treated WGS samples (DR-071).

# Libraries ---------------------------------------------------------------

library(R2CPCT)
library(GenomicRanges)
library(ShatterSeek)

# Initialize a file-logger to store INFO / TRACE messages.
R2CPCT::initializeLogger(file = '~/AbiEnza.log')

# Load metadata of the Abi/Enza-treated patients.
AbiEnza.Metadata <- readxl::read_xlsx('Misc/Suppl. Table 1 - OverviewOfData.xlsx', sheet = 'Sample overview')

# Save metadata.
save(AbiEnza.Metadata, file = '/mnt/onco0002/repository/HMF/DR71/Oct2021/RData/AbiEnza.Metadata.RData')


# Import / convert WGS samples --------------------------------------------

# Import the WGS data of all samples in the cohort.
AbiEnza.CohortWGS <- R2CPCT::importWGSOfCohort(AbiEnza.Metadata$sampleId, '/mnt/onco0002/repository/HMF/DR41/Oct2021/dataHMF/WGS/combinedData/', nThreads = 10)
save(AbiEnza.CohortWGS, file = '/mnt/onco0002/repository/HMF/DR71/Oct2021/RData/AbiEnza.CohortWGS.RData')


# Perform GISTIC2 analysis ------------------------------------------------

# Generate the commands to perform GISTIC2 (perform this in your Bash terminal).
# Perform this separately per responder group (Good / Poor).

# All samples.
R2CPCT::performGISTIC2(AbiEnza.CohortWGS$copyNumbers, outputFolder = '/mnt/onco0002/repository/HMF/DR71/Oct2021/results/GISTIC2/PredictionOfAbiEnza/AllSamples/')

# Good responders.
R2CPCT::performGISTIC2(AbiEnza.CohortWGS$copyNumbers[GenomicRanges::mcols(AbiEnza.CohortWGS$copyNumbers)$sample %in% (AbiEnza.Metadata %>% dplyr::filter(grepl('Good Responder', Responder)) %>% dplyr::pull(sampleId)),], outputFolder = '/mnt/onco0002/repository/HMF/DR71/Oct2021/results/GISTIC2/PredictionOfAbiEnza/GoodResponders/')

# Poor responders.
R2CPCT::performGISTIC2(AbiEnza.CohortWGS$copyNumbers[GenomicRanges::mcols(AbiEnza.CohortWGS$copyNumbers)$sample %in% (AbiEnza.Metadata %>% dplyr::filter(grepl('Bad Responder', Responder)) %>% dplyr::pull(sampleId)),], outputFolder = '/mnt/onco0002/repository/HMF/DR71/Oct2021/results/GISTIC2/PredictionOfAbiEnza/BadResponders/')

# GISTIC2 command used:
#sudo docker run -it --rm -v /mnt/onco0002/repository/HMF/DR71/Oct2021/results/GISTIC2/PredictionOfAbiEnza/:/mnt/onco0002/repository/HMF/DR71/Oct2021/results/GISTIC2/PredictionOfAbiEnza/ -v /tmp/:/tmp/ --entrypoint bash shixiangwang/gistic gistic2 -b /mnt/onco0002/repository/HMF/DR71/Oct2021/results/GISTIC2/PredictionOfAbiEnza/AllSamples/ -seg /tmp/RtmpRZAh8h/60c250100c86.txt -refgene refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat -genegistic 1 -gcm extreme -maxseg 4000 -broad 1 -brlen 0.98 -conf 0.95 -rx 0 -cap 3 -saveseg 0 -armpeel 1 -smallmem 0 -res 0.01 -ta 0.3 -td 0.3 -savedata 0 -savegene 1 -qvt 0.1 -twoside 0
#sudo docker run -it --rm -v /mnt/onco0002/repository/HMF/DR71/Oct2021/results/GISTIC2/PredictionOfAbiEnza/:/mnt/onco0002/repository/HMF/DR71/Oct2021/results/GISTIC2/PredictionOfAbiEnza/ -v /tmp/:/tmp/ --entrypoint bash shixiangwang/gistic gistic2 -b /mnt/onco0002/repository/HMF/DR71/Oct2021/results/GISTIC2/PredictionOfAbiEnza/GoodResponders/ -seg /tmp/RtmpRZAh8h/60c2c1fd0f5.txt -refgene refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat -genegistic 1 -gcm extreme -maxseg 4000 -broad 1 -brlen 0.98 -conf 0.95 -rx 0 -cap 3 -saveseg 0 -armpeel 1 -smallmem 0 -res 0.01 -ta 0.3 -td 0.3 -savedata 0 -savegene 1 -qvt 0.1 -twoside 0
#sudo docker run -it --rm -v /mnt/onco0002/repository/HMF/DR71/Oct2021/results/GISTIC2/PredictionOfAbiEnza/:/mnt/onco0002/repository/HMF/DR71/Oct2021/results/GISTIC2/PredictionOfAbiEnza/ -v /tmp/:/tmp/ --entrypoint bash shixiangwang/gistic gistic2 -b /mnt/onco0002/repository/HMF/DR71/Oct2021/results/GISTIC2/PredictionOfAbiEnza/BadResponders/ -seg /tmp/RtmpRZAh8h/60c2281a0509.txt -refgene refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat -genegistic 1 -gcm extreme -maxseg 4000 -broad 1 -brlen 0.98 -conf 0.95 -rx 0 -cap 3 -saveseg 0 -armpeel 1 -smallmem 0 -res 0.01 -ta 0.3 -td 0.3 -savedata 0 -savegene 1 -qvt 0.1 -twoside 0


# Cohort-wide analysis ----------------------------------------------------

# Initialize a list to contain the results of the WGS analysis.
AbiEnza.Results <- list()

# Generate overview of mutational burden per sample and cohort.
AbiEnza.Results$mutationalBurden <- R2CPCT::determineMutationalBurden(AbiEnza.CohortWGS, minTAF.Muts = 0, minTAF.SV = 0)

# Perform dN/dS on entire cohort and per responder groups.
AbiEnza.Results$dNdS <- R2CPCT::rundNdS(AbiEnza.CohortWGS$somaticVariants)
AbiEnza.Results$dNdS.GoodResponders <- R2CPCT::rundNdS(AbiEnza.CohortWGS$somaticVariants[AbiEnza.CohortWGS$somaticVariants$sample %in% (AbiEnza.Metadata %>% dplyr::filter(grepl('Good Responder', Responder)) %>% dplyr::pull(sampleId)),])
AbiEnza.Results$dNdS.BadResponders <- R2CPCT::rundNdS(AbiEnza.CohortWGS$somaticVariants[AbiEnza.CohortWGS$somaticVariants$sample %in% (AbiEnza.Metadata %>% dplyr::filter(grepl('Bad Responder', Responder)) %>% dplyr::pull(sampleId)),])

# Import the GISTIC2-determined recurrent CNA peaks.
AbiEnza.Results$GISTIC2.AllSamples <- R2CPCT::importGISTIC2('/mnt/onco0002/repository/HMF/DR71/Oct2021/results/GISTIC2/PredictionOfAbiEnza/AllSamples/')
AbiEnza.Results$GISTIC2.GoodResponders <- R2CPCT::importGISTIC2('/mnt/onco0002/repository/HMF/DR71/Oct2021/results/GISTIC2/PredictionOfAbiEnza/GoodResponders/')
AbiEnza.Results$GISTIC2.BadResponders <- R2CPCT::importGISTIC2('/mnt/onco0002/repository/HMF/DR71/Oct2021/results/GISTIC2/PredictionOfAbiEnza/BadResponders/')

# Generate the gene-level mutational overview.
AbiEnza.Results$combinedReport <- R2CPCT::generateCombinedReport(AbiEnza.CohortWGS, dNdS = AbiEnza.Results$dNdS, GISTIC2 = AbiEnza.Results$GISTIC2.AllSamples, nThreads = 20, mutantsOnly = T)

# Determine mutational signatures.
AbiEnza.Results$mutSigs <- R2CPCT::fitMutSigs(AbiEnza.CohortWGS$somaticVariants, restrictiveFit = F)

# Determine # of Ti/Tv and ratio per sample.
AbiEnza.Results$TiTv <- R2CPCT::determineTiTv(AbiEnza.Results$mutSigs$SNV$mutMatrix)

# Save results.
save(AbiEnza.Results, file = '/mnt/onco0002/repository/HMF/DR71/Oct2021/RData/AbiEnza.Results.RData')