# Author:    Job van Riet
# Date:      27-07-22
# Function:  Import and processing of the Abi/Enza-treated RNA-Seq samples (DR-071).

# Libraries and data. ----

library(R2CPCT)
library(DESeq2)

# Load metadata of the Abi/Enza-treated patients.
load('/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.Metadata.RData')

# Load RNA-Seq of the DR-071 cohort.
load('/mnt/share1/repository/HMF/DR71/Dec2021/RData/DR71.RNASeq.RData')

# Retrieve batch-effect genes from the full DR-071 RNA-Seq cohort----
batchGenes <- DR71.RNASeq$DESeq2Results.BetweenMajorBiopsySite %>% dplyr::filter(isSig) %>% dplyr::distinct(ENSEMBL)


# Determine inclusion. ----

# Clean-up the responder category for use in DESeq2.
AbiEnza.Metadata <- AbiEnza.Metadata %>% dplyr::mutate(responderCategory.DESeq2 = gsub(' .*', '', Responder))

# List of training samples.
trainingSamples <- readr::read_delim('Misc/samplesIncludedInTrainingSet.txt', delim = '\t') %>% dplyr::pull(includedTraining)


# Retrieve read-counts of inclusion samples. ----

# Import gene-annotations.
utils::data('GENCODE.v38', package = 'R2CPCT')
geneInfo <- tibble::as_tibble(S4Vectors::mcols(GENCODE.v38))

# Only retain protein-coding genes.
geneInfo <- geneInfo %>% dplyr::filter(gene_type == 'protein_coding')
geneInfo <- geneInfo %>% dplyr::distinct(SYMBOL, ENSEMBL)

# Remove batch-effect genes.
geneInfo <- geneInfo %>% dplyr::filter(!ENSEMBL %in% batchGenes$ENSEMBL)

# Import counts of the entire DR-071 cohort.
countData <- DESeq2::counts(DR71.RNASeq$DESeq2.withoutSequencingBatch, normalized = F, replaced = F)

# Subset genes.
countData <- countData[rownames(countData) %in% geneInfo$ENSEMBL,]

# Subset inclusion samples.
countData.AbiEnza <- countData[,colnames(countData) %in% trainingSamples]


# Generate DESeq2 Dataset. ----

DESeq2Counts.AbiEnza <- DESeq2::DESeqDataSetFromMatrix(countData = countData.AbiEnza, colData = AbiEnza.Metadata[base::match(base::colnames(countData.AbiEnza), AbiEnza.Metadata$hmfSampleId),], design = ~responderCategory.DESeq2)
SummarizedExperiment::rowData(DESeq2Counts.AbiEnza)$ENSEMBL <- BiocGenerics::rownames(DESeq2Counts.AbiEnza)
SummarizedExperiment::rowData(DESeq2Counts.AbiEnza) <- tibble::as_tibble(SummarizedExperiment::rowData(DESeq2Counts.AbiEnza)) %>% dplyr::left_join(geneInfo %>% dplyr::distinct(SYMBOL, ENSEMBL), by = 'ENSEMBL')


# Perform DESeq2. ----

AbiEnza.RNASeq <- list()

AbiEnza.RNASeq$DESeq2 <- DESeq2::DESeq(DESeq2Counts.AbiEnza, test = 'Wald', parallel = F, BPPARAM = BiocParallel::MulticoreParam(workers = 20))

# Retrieve the results. (Bad vs. Good responders)
AbiEnza.RNASeq$DESeq2Results <- R2CPCT::retrieveDESeq2Results(AbiEnza.RNASeq$DESeq2, contrast = c('responderCategory.DESeq2', 'Bad', 'Good'))

# Retrieve Diff. Exprs. genes (DE).
AbiEnza.RNASeq$DESeq2Results <- AbiEnza.RNASeq$DESeq2Results %>% 
  dplyr::mutate(isSig = ifelse(padj <= 0.05 & baseMean >= 25 & lfcSE <= 1 & abs(log2FoldChange) >= 0.5, 'Significant', 'Not Significant'))

save(AbiEnza.RNASeq, file = '/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.RNASeq.RData')
