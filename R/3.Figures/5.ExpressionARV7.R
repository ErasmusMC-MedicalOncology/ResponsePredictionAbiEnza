# Author:    Job van Riet
# Date:      21-03-22
# Function:  Detection of AR-V7.

# Libraries ----

library(R2CPCT)

# Load metadata of the Abi/Enza-treated patients.
load('/mnt/share1/repository/HMF/DR71/Dec2021/RData/AbiEnza.Metadata.RData')

# Retrieve WTS of DR-071.
load('/mnt/share1/repository/HMF/DR71/Dec2021/RData/DR71.RNASeq.RData')

# Visualize and test AR-V7 expression. ----

source('R/3.Figures/misc_Themes.R')

AR.PSI <- DR71.RNASeq$Splicing.AR$PSI %>% 
  dplyr::select(sample, PSI.ARv7, counts.ARv7) %>% 
  dplyr::inner_join(AbiEnza.Metadata, by = c('sample' = 'sampleId'))

# Test.
rstatix::pairwise_wilcox_test(AR.PSI, PSI.ARv7 ~ Responder)

# Plot.
AR.PSI %>% 
  ggplot2::ggplot(., ggplot2::aes(x = Responder, y = PSI.ARv7, fill = Responder)) +
  gghalves::geom_half_boxplot(outlier.shape = NA) +
  gghalves::geom_half_point_panel(shape = 21) +
  ggplot2::scale_fill_manual(values = colorPalette, breaks = unique(AbiEnza.Metadata$Responder), guide = ggplot2::guide_legend(title = NULL, title.position = 'top', title.hjust = 0.5, ncol = 1, keywidth = 0.5, keyheight = 0.5)) +
  ggplot2::scale_y_continuous(breaks = seq(0, 1, by = 0.01), labels = scales::percent_format(accuracy = 1), expand = c(0,0), limits = c(0,.06)) +
  ggplot2::labs(x = 'AR-V7', y = 'AR-V7 PSI (Ψ), relative to ARwt') +
  theme_Job
