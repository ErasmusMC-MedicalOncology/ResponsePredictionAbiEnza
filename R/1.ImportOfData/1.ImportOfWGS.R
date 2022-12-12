# Author:    Job van Riet
# Date:      05-01-22
# Function:  Import and processing of the Abi/Enza-treated WGS samples (DR-071 of Dec. 2021).

# Libraries. ----

library(R2CPCT)
library(GenomicRanges)
library(ShatterSeek)

# Initialize a file-logger to store INFO / TRACE messages.
R2CPCT::initializeLogger(file = '~/AbiEnza.log')

# Import metadata. ----

# Load metadata of the Abi/Enza-treated patients.
AbiEnza.Metadata <- readxl::read_xlsx('Misc/Suppl. Table 1 - OverviewOfData.xlsx', sheet = 'Sample Information')

# Save metadata.
save(AbiEnza.Metadata, file = '/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.Metadata.RData')


# Import WGS data with R2CPCT. ----

# Import the WGS data of all samples in the cohort.
AbiEnza.CohortWGS <- R2CPCT::importWGSOfCohort(AbiEnza.Metadata$sampleId, '/mnt/share1/repository/HMF/DR71/Dec2021/dataHMF/combinedData/', nThreads = 8)
save(AbiEnza.CohortWGS, file = '/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.CohortWGS.RData')


# Perform GISTIC2. -----

# Generate the commands to perform GISTIC2 (perform this in your Bash terminal with the GISTIC2 docker).
# Perform this separately per responder category (Good / Bad).

# All samples.
R2CPCT::performGISTIC2(AbiEnza.CohortWGS$copyNumbers, outputFolder = '/mnt/share1/repository/HMF/DR71/Dec2021/results/GISTIC2/PredictionOfAbiEnza/AllSamples/')

# Good responders.
R2CPCT::performGISTIC2(AbiEnza.CohortWGS$copyNumbers[GenomicRanges::mcols(AbiEnza.CohortWGS$copyNumbers)$sample %in% (AbiEnza.Metadata %>% dplyr::filter(grepl('Good Responder', Responder)) %>% dplyr::pull(sampleId)),], outputFolder = '/mnt/share1/repository/HMF/DR71/Dec2021/results/GISTIC2/PredictionOfAbiEnza/GoodResponders/')

# Poor responders.
R2CPCT::performGISTIC2(AbiEnza.CohortWGS$copyNumbers[GenomicRanges::mcols(AbiEnza.CohortWGS$copyNumbers)$sample %in% (AbiEnza.Metadata %>% dplyr::filter(grepl('Bad Responder', Responder)) %>% dplyr::pull(sampleId)),], outputFolder = '/mnt/share1/repository/HMF/DR71/Dec2021/results/GISTIC2/PredictionOfAbiEnza/BadResponders/')

# GISTIC2 command used:
# Repeated samples in separate analysis

#sudo docker run -it --rm -v /mnt/share1/repository/HMF/DR71/Dec2021/results/GISTIC2/PredictionOfAbiEnza/:/mnt/share1/repository/HMF/DR71/Dec2021/results/GISTIC2/PredictionOfAbiEnza/ -v /tmp/:/tmp/ --entrypoint bash shixiangwang/gistic gistic2 -b /mnt/share1/repository/HMF/DR71/Dec2021/results/GISTIC2/PredictionOfAbiEnza/AllSamples/ -seg /tmp/RtmpyYjjN7/3ae358130a68.txt -refgene refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat -genegistic 1 -gcm extreme -maxseg 4000 -broad 1 -brlen 0.98 -conf 0.95 -rx 0 -cap 3 -saveseg 0 -armpeel 1 -smallmem 0 -res 0.01 -ta 0.3 -td 0.3 -savedata 0 -savegene 1 -qvt 0.1 -twoside 0
#sudo docker run -it --rm -v /mnt/share1/repository/HMF/DR71/Dec2021/results/GISTIC2/PredictionOfAbiEnza/:/mnt/share1/repository/HMF/DR71/Dec2021/results/GISTIC2/PredictionOfAbiEnza/ -v /tmp/:/tmp/ --entrypoint bash shixiangwang/gistic gistic2 -b /mnt/share1/repository/HMF/DR71/Dec2021/results/GISTIC2/PredictionOfAbiEnza/GoodResponders/ -seg /tmp/RtmpaJ7Bk5/4eea3d339bc9.txt -refgene refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat -genegistic 1 -gcm extreme -maxseg 4000 -broad 1 -brlen 0.98 -conf 0.95 -rx 0 -cap 3 -saveseg 0 -armpeel 1 -smallmem 0 -res 0.01 -ta 0.3 -td 0.3 -savedata 0 -savegene 1 -qvt 0.1 -twoside 0
#sudo docker run -it --rm -v /mnt/share1/repository/HMF/DR71/Dec2021/results/GISTIC2/PredictionOfAbiEnza/:/mnt/share1/repository/HMF/DR71/Dec2021/results/GISTIC2/PredictionOfAbiEnza/ -v /tmp/:/tmp/ --entrypoint bash shixiangwang/gistic gistic2 -b /mnt/share1/repository/HMF/DR71/Dec2021/results/GISTIC2/PredictionOfAbiEnza/BadResponders/ -seg /tmp/RtmpaJ7Bk5/4eeac42c45c.txt -refgene refgenefiles/hg19.UCSC.add_miR.140312.refgene.mat -genegistic 1 -gcm extreme -maxseg 4000 -broad 1 -brlen 0.98 -conf 0.95 -rx 0 -cap 3 -saveseg 0 -armpeel 1 -smallmem 0 -res 0.01 -ta 0.3 -td 0.3 -savedata 0 -savegene 1 -qvt 0.1 -twoside 0


# Cohort-wide analysis ----------------------------------------------------

# Initialize a list to contain the results of the WGS analysis.
AbiEnza.Results <- list()

# Generate overview of mutational burden per sample and cohort.
AbiEnza.Results$mutationalBurden <- R2CPCT::determineMutationalBurden(AbiEnza.CohortWGS, minTAF.Muts = 0, minTAF.SV = 0)

# Determine exonic TMB.
data(GENCODE.v38, package = 'R2CPCT')
GENCODE.v38.proteinCoding <- GENCODE.v38[GENCODE.v38$gene_type == 'protein_coding']
GENCODE.v38.proteinCoding <- GenomicRanges::reduce(GENCODE.v38.proteinCoding)

exonMuts <- IRanges::subsetByOverlaps(AbiEnza.CohortWGS$somaticVariants, GENCODE.v38.proteinCoding) %>% 
    tibble::as_tibble(S4Vectors::mcols(.)) %>% 
    dplyr::group_by(sample) %>% 
    dplyr::summarise(
        Exome.TMB = dplyr::n() / (sum(GenomicRanges::width(GENCODE.v38.proteinCoding))/1e+06)
    ) %>% 
    dplyr::ungroup()

AbiEnza.Results$mutationalBurden <- AbiEnza.Results$mutationalBurden %>% dplyr::inner_join(exonMuts)

# Perform dN/dS on entire cohort and per responder groups.
AbiEnza.Results$dNdS <- R2CPCT::rundNdS(AbiEnza.CohortWGS$somaticVariants)
AbiEnza.Results$dNdS.GoodResponders <- R2CPCT::rundNdS(AbiEnza.CohortWGS$somaticVariants[AbiEnza.CohortWGS$somaticVariants$sample %in% (AbiEnza.Metadata %>% dplyr::filter(grepl('Good Responder', Responder)) %>% dplyr::pull(sampleId)),])
AbiEnza.Results$dNdS.BadResponders <- R2CPCT::rundNdS(AbiEnza.CohortWGS$somaticVariants[AbiEnza.CohortWGS$somaticVariants$sample %in% (AbiEnza.Metadata %>% dplyr::filter(grepl('Bad Responder', Responder)) %>% dplyr::pull(sampleId)),])

# Import the GISTIC2-determined recurrent CNA peaks.
AbiEnza.Results$GISTIC2.AllSamples <- R2CPCT::importGISTIC2('/mnt/share1/repository/HMF/DR71/Dec2021/results/GISTIC2/PredictionOfAbiEnza/AllSamples/')
AbiEnza.Results$GISTIC2.GoodResponders <- R2CPCT::importGISTIC2('/mnt/share1/repository/HMF/DR71/Dec2021/results/GISTIC2/PredictionOfAbiEnza/GoodResponders/')
AbiEnza.Results$GISTIC2.BadResponders <- R2CPCT::importGISTIC2('/mnt/share1/repository/HMF/DR71/Dec2021/results/GISTIC2/PredictionOfAbiEnza/BadResponders/')

# Generate the gene-level mutational overview.
AbiEnza.Results$combinedReport <- R2CPCT::generateCombinedReport(AbiEnza.CohortWGS, dNdS = AbiEnza.Results$dNdS, GISTIC2 = AbiEnza.Results$GISTIC2.AllSamples, nThreads = 40, mutantsOnly = T)

# Determine mutational signatures.
AbiEnza.Results$mutSigs <- R2CPCT::fitMutSigs(AbiEnza.CohortWGS$somaticVariants, restrictiveFit = T)

# Determine # of Ti/Tv and ratio per sample.
AbiEnza.Results$TiTv <- R2CPCT::determineTiTv(AbiEnza.Results$mutSigs$SNV$mutMatrix)

# Save results.
save(AbiEnza.Results, file = '/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.Results.RData')