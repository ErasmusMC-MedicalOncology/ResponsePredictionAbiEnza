# All utilized and custom code used in the Response Predictor to Abi/Enza.

This repository contains the custom code (R) used to perform all analysis of the study as further detailed by De Jong and Danyi et al. in <JOURNAL> (2021): [TITLE](https://www.google.com/).

All processed (and cleaned) data as-used in the presented analysis and figures has been published alongside the manuscript (Suppl. Table 1).

The raw and pre-processed data of the CPCT-02 WGS and RNA-Seq cohort can be requested under the data-request **DR-071** from the Hartwig Medical Foundation (HMF): https://www.hartwigmedicalfoundation.nl/en/applying-for-data/. The WGS and RNA-Seq was further annotated and post-processed using the [R2CPCT](https://github.com/J0bbie/R2CPCT) package (v0.4).

The commands and tools used in analysing the RNA-Seq data (starting from paired-end .fastq) are described within the manuscript. This RNA-Seq workflow can be found in ...

All provided code in the R/ folder was performed on Ubuntu 18.04.5 LTS with the following *R* version:
```R
R version 4.0.4 (2021-02-15) -- "Lost Library Book"
Copyright (C) 2021 The R Foundation for Statistical Computing
Platform: x86_64-pc-linux-gnu (64-bit)
```

Required R dependencies are listed on the top of each script.