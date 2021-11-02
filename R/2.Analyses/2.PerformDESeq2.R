# Author:    Job van Riet
# Date:      04-06-21
# Function:  Performs diff. exprs. analysis between poor and good responders.

# Libraries ---------------------------------------------------------------

library(R2CPCT)
library(DESeq2)

# Retrieve RNA-Seq counts.
load('/mnt/data2/hartwig/DR71/Apr2021_AbiEnza/RData/DESeq2Counts.AbiEnza.RData')


# Perform DESeq2 ----------------------------------------------------------

DESeq2.AbiEnza <- DESeq2::DESeq(DESeq2Counts.AbiEnza, test = 'Wald', parallel = T, BPPARAM = BiocParallel::MulticoreParam(workers = 20))

# Retrieve the results. (Poor vs. Good responders)
DESeq2.AbiEnza.Results <- R2CPCT::retrieveDESeq2Results(DESeq2.AbiEnza, contrast = c('responderCategory.DESeq2', 'Bad', 'Good'))

# Retrieve Diff. Exprs. genes (DE).
DESeq2.AbiEnza.Results <- DESeq2.AbiEnza.Results %>% 
  dplyr::mutate(isSig = ifelse(padj <= 0.05 & baseMean >= 25 & lfcSE <= 1 & abs(log2FoldChange) >= 0.5, 'Significant', 'Not Significant'))

# Write to tmp. file and add to Suppl. Table 1.
write.table(DESeq2.AbiEnza.Results, file = 'DifferentialAnalysis.txt', quote = F, sep = '\t', row.names = F)
