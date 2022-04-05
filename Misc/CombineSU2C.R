meta <- readr::read_csv('~/Downloads/su2c/su2c_444_metadata.csv')
genes <- readr::read_csv('~/Downloads/su2c/su2c_444_genes.csv')

counts <- data.frame(readr::read_csv('~/Downloads/su2c/su2c_444_counts.csv'), check.names = F)
rownames(counts) <- genes$gene_id

fpkm <- data.frame(readr::read_csv('~/Downloads/su2c/su2c_444_fpkm.csv'), check.names = F)
rownames(fpkm) <- sprintf('%s_%s', genes$gene_name, genes$gene_id)

# Combine into a list.
data.SU2C <- list()

# Add raw counts/metadata into DESeq2 object.
data.SU2C$Counts <- DESeq2::DESeqDataSetFromMatrix(countData = counts, colData = meta, design = ~assay)
SummarizedExperiment::rowData(data.SU2C$Counts) <- genes

# Add FPKM
data.SU2C$FPKM <- fpkm

# Save RData object.
save(data.SU2C, file = '/mnt/share1/repository/HMF/DR71/Misc/SU2C/SU2C.RData')