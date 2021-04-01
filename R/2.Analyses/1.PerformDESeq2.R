# Author:    Job van Riet
# Date:      01-04-21
# Function:  Performs diff. exprs. analysis between poor and good responders.

# Libraries ---------------------------------------------------------------

library(R2CPCT)
library(DESeq2)

# Load metadata of the Abi/Enza-treated patients.
AbiEnza.Metadata <- readxl::read_xlsx('Misc/Suppl. Table 1 - OverviewOfData.xlsx', sheet = 'Sample overview')

# Retrieve RNA-Seq counts.
load('/mnt/data2/hartwig/DR71/Oct2020_AbiEnza/RData/DESeq2Counts.AbiEnza.RData')

# Perform DESeq2 ----------------------------------------------------------

DESeq2.AbiEnza <- DESeq2::DESeq(DESeq2Counts.AbiEnza, test = 'Wald', parallel = T, BPPARAM = BiocParallel::MulticoreParam(workers = 20))

# Retrieve the results. (Poor vs. Good responders)
DESeq2.AbiEnza.Results <- R2CPCT::retrieveDESeq2Results(DESeq2.AbiEnza, contrast = c('responderCategory.DESeq2', 'Poor', 'Good'))

# Retrieve Diff. Exprs. genes (DE).
DESeq2.AbiEnza.Results <- DESeq2.AbiEnza.Results %>% 
  dplyr::mutate(isSig = ifelse(padj <= 0.05 & lfcSE <= 1 & baseMean >= 50 & abs(log2FoldChange) >= 0.5, 'Significant', 'Not Significant'))

# Write to tmp. file and add to Suppl. Table 1.
write.table(DESeq2.AbiEnza.Results, file = 'DifferentialAnalysis.txt', quote = F, sep = '\t', row.names = F)
