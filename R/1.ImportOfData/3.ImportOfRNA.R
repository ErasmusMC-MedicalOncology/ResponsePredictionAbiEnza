# Author:    Job van Riet
# Date:      26-05-21
# Function:  Import and processing of the Abi/Enza-treated RNA-Seq samples (DR-071).


# Libraries ---------------------------------------------------------------

library(R2CPCT)
library(DESeq2)

# Load metadata of the Abi/Enza-treated patients.
load('/mnt/data2/hartwig/DR71/Apr2021_AbiEnza/RData/AbiEnza.Metadata.RData')

# Clean-up the responder category for use in DESEq2.
AbiEnza.Metadata <- AbiEnza.Metadata %>% dplyr::mutate(responderCategory.DESeq2 = gsub(' .*', '', Responder))

# Select duplicate RNA-Seq samples.
load('/mnt/data2/hartwig/DR71/Apr2021/RData/DR71.MetaData.RData')
repeatedBiopsies <- DR71.MetaData$sampleInfo %>% 
  dplyr::filter(hasMatchingRNA == 'Yes', hmfPatientId %in% AbiEnza.Metadata$hmfPatientId) %>% 
  dplyr::group_by(hmfPatientId) %>% 
  dplyr::mutate(totalRNA = dplyr::n_distinct(hmfSampleId)) %>% 
  dplyr::filter(totalRNA > 1) %>% 
  dplyr::ungroup()


# Retrieve batch-effect genes ---------------------------------------------

load('/mnt/data2/hartwig/DR71/Oct2020/RData/DR71.RNASeq.RData')

## Retrieve batch-effect genes ----
diffGenes.Batch <- DR71.RNASeq$DESeq2Results.BetweenMajorBiopsySite %>% dplyr::filter(padj <= 0.05, lfcSE <= 1, abs(log2FoldChange) >= 1) %>% dplyr::distinct(ENSEMBL)


# Import read-counts ------------------------------------------------------

# Import gene-annotations.
geneInfo <- R2CPCT::GENCODE.v35
geneInfo <- tibble::as_tibble(S4Vectors::mcols(geneInfo))

# Only retain protein-coding genes.
geneInfo <- geneInfo %>% dplyr::filter(gene_type == 'protein_coding')
geneInfo <- geneInfo %>% dplyr::distinct(SYMBOL, ENSEMBL)

# Remove batch-effect genes.
geneInfo <- geneInfo %>% dplyr::filter(!ENSEMBL %in% diffGenes.Batch$ENSEMBL)

# Import counts of the entire DR-071 cohort.
counts <- readr::read_delim('/mnt/data2/hartwig/DR71/Apr2021/dataHMF/RNASeq/counts/DR71_RNA.counts', delim = '\t', comment = '#')
base::colnames(counts) <- base::gsub('_.*', '', base::gsub('.*/', '', base::colnames(counts)))

# Convert to matrix.
countMatrix <- as.matrix(counts[7:ncol(counts)])
rownames(countMatrix) <- counts$Geneid

# Subset Abi/Enza and repeated biopsies samples.
countMatrix <- countMatrix[,colnames(countMatrix) %in% c(repeatedBiopsies$sample, AbiEnza.Metadata$sampleId)]

# Subset genes.
countMatrix <- countMatrix[rownames(countMatrix) %in% geneInfo$ENSEMBL,]


# Generate DESeq2 Dataset -------------------------------------------------

cleanColData <- AbiEnza.Metadata[match(colnames(countMatrix), AbiEnza.Metadata$sampleId),]
cleanColData <- cleanColData %>% dplyr::mutate(responderCategory.DESeq2 = ifelse(is.na(responderCategory.DESeq2), 'Repeat', responderCategory.DESeq2))

DESeq2Counts.AbiEnza <- DESeq2::DESeqDataSetFromMatrix(countData = countMatrix, colData = cleanColData, design = ~responderCategory.DESeq2)
SummarizedExperiment::rowData(DESeq2Counts.AbiEnza)$ENSEMBL <- BiocGenerics::rownames(DESeq2Counts.AbiEnza)
SummarizedExperiment::rowData(DESeq2Counts.AbiEnza) <- tibble::as_tibble(SummarizedExperiment::rowData(DESeq2Counts.AbiEnza)) %>% dplyr::left_join(tibble::as_tibble(S4Vectors::mcols(R2CPCT::GENCODE.v35)) %>% dplyr::distinct(SYMBOL, ENSEMBL), by = 'ENSEMBL')

save(DESeq2Counts.AbiEnza, file = '/mnt/data2/hartwig/DR71/Apr2021_AbiEnza/RData/DESeq2Counts.AbiEnza_NoRepeatedBiopsies.RData')
