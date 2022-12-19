# Author:                       Job van Riet
# Date of modification:         31-01-22
# Function:                     Perform GSEA between Good vs. Poor responders.


# Load packages -----------------------------------------------------------

library(plyr)
library(dplyr)
library(DESeq2)

# Import data -------------------------------------------------------------

# Load metadata of the Abi/Enza-treated patients.
load('/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.Metadata.RData')

# Retrieve WTS-data.
load('/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.RNASeq.RData')


# Define gene-sets --------------------------------------------------------

## Retrieve gene-sets. ----
genesets <- list()

genesets$HALLMARK <- hypeR::msigdb_gsets("Homo sapiens", "H", "", clean = T)$genesets
names(genesets$HALLMARK) <- paste0(names(genesets$HALLMARK), ' (Hallmark)')

genesets$WIKI <- hypeR::msigdb_gsets(species='Homo sapiens', category='C2', subcategory='CP:WIKIPATHWAYS', clean = T)$genesets
names(genesets$WIKI) <- paste0(names(genesets$WIKI), ' (WikiPathways)')

genesets$REACTOME <- hypeR::msigdb_gsets(species='Homo sapiens', category='C2', subcategory='CP:REACTOME', clean = T)$genesets
names(genesets$REACTOME) <- paste0(names(genesets$REACTOME), ' (Reactome)')

genesets$KEGG <- hypeR::msigdb_gsets(species='Homo sapiens', category='C2', subcategory='CP:KEGG', clean = T)$genesets
names(genesets$KEGG) <- paste0(names(genesets$KEGG), ' (KEGG)')

# Select which gene-sets to look for.
testSets <- c(genesets$HALLMARK, genesets$WIKI)
testSets <- testSets[S4Vectors::elementNROWS(testSets) >= 10]
testSets <- testSets[S4Vectors::elementNROWS(testSets) <= 300]

testSets <- testSets[!duplicated(gsub(' \\(.*', '', names(testSets)))]


# Perform GSEA ------------------------------------------------------------

performGSEA <- function(data){
  
  # Retrieve t-statistics.
  y <- data %>% dplyr::filter(!is.na(SYMBOL), !is.na(stat)) %>% 
    dplyr::filter(baseMean >= 10) %>% 
    dplyr::distinct() %>% 
    dplyr::group_by(SYMBOL) %>% 
    dplyr::summarize(stat=mean(stat))
  
  # Sort on t.
  y.ordered <- (y %>% dplyr::select(SYMBOL, stat) %>% dplyr::arrange(stat)) %>% tibble::deframe()
  
  # Run fgsea
  fgseaRes <- fgsea::fgseaMultilevel(
    pathways = testSets, 
    stats    = y.ordered,
    minSize  = 10,
    maxSize  = 300,
    eps = 0
  )
  
  # Add the contrast to keep track of input origin.
  fgseaRes$contrast <- unique(data$contrast)
  
  # Return.
  return(fgseaRes)
  
}

AbiEnza.RNASeq$GSEA <- performGSEA(AbiEnza.RNASeq$DESeq2Results)

# Save --------------------------------------------------------------------

save(AbiEnza.RNASeq, file = '/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.RNASeq.RData')
