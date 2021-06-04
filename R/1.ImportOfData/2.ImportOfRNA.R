# Author:    Job van Riet
# Date:      04-06-21
# Function:  Import and processing of the Abi/Enza-treated RNA-Seq samples (DR-071).


# Libraries ---------------------------------------------------------------

library(R2CPCT)
library(DESeq2)

# Load metadata of the Abi/Enza-treated patients.
load('/mnt/data2/hartwig/DR71/Apr2021_AbiEnza/RData/AbiEnza.Metadata.RData')

# Load metadata of the DR-071 cohort.
load('/mnt/data2/hartwig/DR71/Apr2021/RData/DR71.MetaData.RData')

# Load RNA-Seq of the DR-071 cohort.
load('/mnt/data2/hartwig/DR71/Apr2021/RData/DR71.RNASeq.RData')


# Determine samples -------------------------------------------------------

# Clean-up the responder category for use in DESEq2.
AbiEnza.Metadata <- AbiEnza.Metadata %>% dplyr::mutate(responderCategory.DESeq2 = gsub(' .*', '', Responder))

# Additional repeated biopsies.
repeatedBiopsies <- DR71.MetaData$sampleInfo %>% 
  dplyr::filter(sample %in% c('DRUP01050020T', 'DRUP01070071T', 'DRUP01070057T', 'DRUP01010097T', 'CPCT02070055TII', 'CPCT02020351TII', 'CPCT02070107TII', 'CPCT02010692TII', 'CPCT02140041TII', 'CPCT02140041TIII')) %>% 
  dplyr::select(sample, hmfSampleId)


# Retrieve batch-effect genes ---------------------------------------------

## Retrieve batch-effect genes ----
batchGenes <- DR71.RNASeq$DESeq2Results.BetweenMajorBiopsySite %>% dplyr::filter(isSig) %>% dplyr::distinct(ENSEMBL)


# Retrieve read-counts ----------------------------------------------------

# Import gene-annotations.
geneInfo <- R2CPCT::GENCODE.v35
geneInfo <- tibble::as_tibble(S4Vectors::mcols(geneInfo))

# Only retain protein-coding genes.
geneInfo <- geneInfo %>% dplyr::filter(gene_type == 'protein_coding')
geneInfo <- geneInfo %>% dplyr::distinct(SYMBOL, ENSEMBL)

# Remove batch-effect genes.
geneInfo <- geneInfo %>% dplyr::filter(!ENSEMBL %in% batchGenes$ENSEMBL)

# Import counts of the entire DR-071 cohort.
countData <- DESeq2::counts(DR71.RNASeq$DESeq2.withoutSequencingBatch, normalized = F, replaced = F)

# Subset Abi/Enza and repeated biopsies samples.
countData <- countData[,colnames(countData) %in% c(repeatedBiopsies$hmfSampleId, AbiEnza.Metadata$hmfSampleId)]

# Subset genes.
countData <- countData[rownames(countData) %in% geneInfo$ENSEMBL,]


# Generate DESeq2 Dataset -------------------------------------------------

cleanColData <- AbiEnza.Metadata[match(colnames(countData), AbiEnza.Metadata$hmfSampleId),]
cleanColData <- cleanColData %>% dplyr::mutate(responderCategory.DESeq2 = ifelse(is.na(responderCategory.DESeq2), 'Repeat', responderCategory.DESeq2))

DESeq2Counts.AbiEnza <- DESeq2::DESeqDataSetFromMatrix(countData = countData, colData = cleanColData, design = ~responderCategory.DESeq2)
SummarizedExperiment::rowData(DESeq2Counts.AbiEnza)$ENSEMBL <- BiocGenerics::rownames(DESeq2Counts.AbiEnza)
SummarizedExperiment::rowData(DESeq2Counts.AbiEnza) <- tibble::as_tibble(SummarizedExperiment::rowData(DESeq2Counts.AbiEnza)) %>% dplyr::left_join(tibble::as_tibble(S4Vectors::mcols(R2CPCT::GENCODE.v35)) %>% dplyr::distinct(SYMBOL, ENSEMBL), by = 'ENSEMBL')

save(DESeq2Counts.AbiEnza, file = '/mnt/data2/hartwig/DR71/Apr2021_AbiEnza/RData/DESeq2Counts.AbiEnza.RData')
