# Author:    Job van Riet
# Date:      26-01-23
# Function:  Import and processing of the Abi/Enza-treated RNA-Seq samples (DR-071).

# Libraries and data. ----

library(DESeq2)
library(dplyr)

# Load WTS of the DR-071 cohort.
load('/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.RNASeq.RData')


# Perform DE in LOOCV fashion. ----

# Loop over each sample and remove one sample prior to Poor vs. Good DE to perform DE in a LOOCV fashion.
DE.LOOCV <- pbapply::pblapply(AbiEnza.RNASeq$DESeq2$hmfSampleId, function(sampleOut){
    
    x <- DESeq2::DESeq(AbiEnza.RNASeq$DESeq2[,AbiEnza.RNASeq$DESeq2$hmfSampleId != sampleOut], test = 'Wald', parallel = F, BPPARAM = BiocParallel::MulticoreParam(workers = 8))
    
    # Get diff. results.
    results <- DESeq2::results(x, pAdjustMethod = 'BH', filterFun = IHW::ihw, parallel = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 8), tidy = FALSE, contrast = c('responderCategory.DESeq2', 'Bad', 'Good'))
    
    # Shrink the LFC.
    results.LFC <- DESeq2::lfcShrink(x, res = results, parallel = TRUE, BPPARAM = BiocParallel::MulticoreParam(workers = 8), contrast = c('responderCategory.DESeq2', 'Bad', 'Good'), type = 'ashr')
    
    # Add ENSEMBL as column.
    results.LFC$ENSEMBL <- BiocGenerics::rownames(results.LFC)
    
    # Re-add the t-statistic
    results.LFC$stat <- results[match(rownames(results), rownames(results.LFC)),]$stat
    
    # Re-add the q weight.
    results.LFC$weight <- results[match(rownames(results), rownames(results.LFC)),]$weight
    
    # Convert to tibble.
    results.LFC <- tibble::as_tibble(results.LFC)
    
    # Add contrast.
    results.LFC$contrast <- 'Poor vs. Good Responders'
    
    # Add the SYMBOL.
    results.LFC <- results.LFC %>% dplyr::left_join(tibble::as_tibble(rowData(x)[c('ENSEMBL', 'SYMBOL')]) %>% dplyr::distinct(SYMBOL, ENSEMBL), by = 'ENSEMBL')
    
    results.LFC$sampleOut <- factor(sampleOut)
    
    results.LFC <- results.LFC %>% dplyr::mutate(
        isSig = ifelse(padj <= 0.05 & baseMean >= 25 & lfcSE <= 1 & abs(log2FoldChange) >= 0.5, 'Significant', 'Not Significant'),
        sampleOut = factor(sampleOut)
    )
    
    return(results.LFC)
    
}, cl = 10)

# Determine recurrent DEGs.
DE.LOOCV <- dplyr::bind_rows(DE.LOOCV) %>% dplyr::mutate_if(is.character,as.factor)
DE.LOOCV.Summed <- DE.LOOCV %>% dplyr::filter(isSig == 'Significant') %>% dplyr::group_by(ENSEMBL, SYMBOL) %>% dplyr::summarise(totalSig = dplyr::n_distinct(sampleOut)) %>% dplyr::ungroup()

# Save to object.
DE.LOOCV <- list(DE.LOOCV = DE.LOOCV, DE.LOOCV.Summed = DE.LOOCV.Summed)
saveRDS(DE.LOOCV, file = '/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.LOOCV.Rds')
