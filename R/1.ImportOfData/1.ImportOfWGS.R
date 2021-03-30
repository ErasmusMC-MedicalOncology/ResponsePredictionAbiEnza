# Author:    Job van Riet
# Date:      30-03-21
# Function:  Import and processing of the Abi/Enza-treated WGS samples (DR-071).

# Libraries ---------------------------------------------------------------

library(R2CPCT)
library(GenomicRanges)
library(ShatterSeek)

# Initialize a file-logger to store INFO / TRACE messages.
R2CPCT::initializeLogger(file = '~/AbiEnza.log')

# Load metadata of the Abi/Enza-treated patients.
AbiEnza.Metadata <- readxl::read_xlsx('Misc/SupTable1_OverviewOfData.xlsx', sheet = 'Sample overview')


# Retrieve all relevant HMF files -----------------------------------------

