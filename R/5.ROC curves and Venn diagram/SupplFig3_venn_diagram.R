library(ggplot2)
library(RColorBrewer)
library(patchwork)
require(gridExtra)
library(grid)
library(ggvenn)
library(readxl)
library(caret)
library(plyr)


# data <- read_excel("./Predictions_HMF.xlsx", sheet = "Predictions - HMF")
# data <- data %>% dplyr::filter(subgroupCohort == 'Training')
# data_df <- as.data.frame(data)
#
# sample_info <- read_excel("./Suppl.\ Table\ 1\ -\ OverviewOfData.xlsx", sheet = "Sample information")
# sample_info <- sample_info[sample_info$hmfSampleId %in% data_df$hmfSampleId, ]
# sample_info$Responder[sample_info$Responder == "Bad Responder (≤100 days)"] <- 0
# sample_info$Responder[sample_info$Responder == "Good Responder (≥180 days)"] <- 1
#
# pred_data <- data_df[, c('hmfSampleId', 'subgroupCohort', 'WGS_probability_Poor', 'WGS_probability_Good', 'Transcriptomics_probability_Poor', 'Transcriptomics_probability_Good')]
#
# pred_data <- pred_data[order(pred_data$hmfSampleId),]
# sample_info <- sample_info[order(sample_info$hmfSampleId),]
#
# responder_ground_truth <- factor(sample_info$Responder)
# pred_data$Transcriptomics_probability_Good[pred_data$Transcriptomics_probability_Good > 0.5] <- 1
# pred_data$Transcriptomics_probability_Good[pred_data$Transcriptomics_probability_Good <= 0.5] <- 0
# pred_data$WGS_probability_Good[pred_data$WGS_probability_Good > 0.5] <- 1
# pred_data$WGS_probability_Good[pred_data$WGS_probability_Good <= 0.5] <- 0
#
# transcriptomics_pred <- factor(pred_data$Transcriptomics_probability_Good)
# genomics_pred <- factor(pred_data$WGS_probability_Good)
#
# transcriptomics_confusionmatrix <- confusionMatrix(responder_ground_truth, transcriptomics_pred)
# genomics_confusionmatrix <- confusionMatrix(responder_ground_truth, genomics_pred)


#### TP (true positives), FP (false positives), TN (true negatives) and FN (false negatives) can be extracted from
#### the Suppl. Table 1 / Source Data
#### predictions for WGS-only and WTS-only models: "Predictions - HMF" sheet
#### ground truth labels: "Sample information" sheet ("Responder")

png(file="./Fig_venn_diagram.png", width=9, height = 9, units = "in", res= 900, pointsize = 4)


data <- read.csv("./Venn_predictions_tp.csv", header = TRUE, sep = ",")
data$genomic = as.logical(data$genomic)
data$transcriptomic = as.logical(data$transcriptomic)
# ,sample_ids,genomic,transcriptomic
# 0,HMF001410A,1,0
# 1,HMF002337A,1,1
# 2,HMF000879A,1,1
# 3,HMF002964A,1,1
# 4,HMF000574A,0,1
# 5,HMF003830A,1,1
# 6,HMF000679A,0,1
# 7,HMF002427A,1,1
# 8,HMF000050A,1,0
# 9,HMF001189A,1,0
# 10,HMF000599A,1,1
# 11,HMF004001A,0,1
# 12,HMF004373A,1,1
# 13,HMF004221A,1,0
# 14,HMF001344B,1,1
# 15,HMF004443A,1,1
# 16,HMF004355A,1,0
# 17,HMF001378B,1,0
# 18,HMF000377A,1,0
# 19,HMF002145A,1,1
# 20,HMF002564A,1,1
# 21,HMF000652A,1,1
# 22,HMF000822A,1,0
# 23,HMF002674A,1,1
# 24,HMF002476A,1,1
# 25,HMF003789A,1,1
# 26,HMF000986A,1,0
# 27,HMF003826A,1,1
# 28,HMF004057A,0,1
# 29,HMF004363A,1,0
# 30,HMF003914A,1,1
# 31,HMF004387A,1,0
# 32,HMF003728A,1,1
# 33,HMF000759A,1,0


plot1 <- ggvenn(data,
  fill_color = c("#0073C2FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4, show_percentage=FALSE
  ) + labs(title= "True Positives") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
      strip.background = ggplot2::element_rect(colour = 'white', fill = 'white'),
      panel.background = element_rect(fill='white', colour='white'),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank())

data_fp <- read.csv("./Venn_predictions_fp.csv", header = TRUE, sep = ",")
data_fp$genomic = as.logical(data_fp$genomic)
data_fp$transcriptomic = as.logical(data_fp$transcriptomic)
# ,sample_ids,genomic,transcriptomic
# 0,HMF000609A,1,0
# 1,HMF003638A,1,1
# 2,HMF003192A,0,1
# 3,HMF001412A,1,0
# 4,HMF004113A,1,0
# 5,HMF000376B,1,0
# 6,HMF004216A,1,0
# 7,HMF001186A,0,1
# 8,HMF002279A,1,0
# 9,HMF001758A,1,0
# 10,HMF002177A,1,1
# 11,HMF002580A,1,0
# 12,HMF004351A,1,0
# 13,HMF003918A,1,0
# 14,HMF004109A,1,0
# 15,HMF003628A,1,0
# 16,HMF000938A,1,1
# 17,HMF001909A,1,0
# 18,HMF004083A,1,0
# 19,HMF003821A,0,1
# 20,HMF002026A,1,0
# 21,HMF002672A,1,0
# 22,HMF004234A,1,1
# 23,HMF003990A,1,0

plot2 <- ggvenn(data_fp,
  fill_color = c("#0073C2FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4, show_percentage=FALSE
  ) + labs(title= "False Positives") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
      strip.background = ggplot2::element_rect(colour = 'white', fill = 'white'),
      panel.background = element_rect(fill='white', colour='white'),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank())

data_tn <- read.csv("./Venn_predictions_tn.csv", header = TRUE, sep = ",")
data_tn$genomic = as.logical(data_tn$genomic)
data_tn$transcriptomic = as.logical(data_tn$transcriptomic)
# ,sample_ids,genomic,transcriptomic
# 0,HMF000332A,1,1
# 1,HMF001024A,1,1
# 2,HMF001885A,1,1
# 3,HMF000609A,0,1
# 4,HMF000983A,1,1
# 5,HMF001697A,1,1
# 6,HMF001925B,1,1
# 7,HMF000384A,1,1
# 8,HMF003192A,1,0
# 9,HMF001412A,0,1
# 10,HMF004113A,0,1
# 11,HMF003899A,1,1
# 12,HMF000376B,0,1
# 13,HMF004216A,0,1
# 14,HMF002207A,1,1
# 15,HMF000222A,1,1
# 16,HMF001186A,1,0
# 17,HMF002279A,0,1
# 18,HMF001758A,0,1
# 19,HMF002580A,0,1
# 20,HMF004351A,0,1
# 21,HMF003918A,0,1
# 22,HMF004189A,1,1
# 23,HMF004109A,0,1
# 24,HMF001337A,1,1
# 25,HMF003628A,0,1
# 26,HMF004348A,1,1
# 27,HMF001909A,0,1
# 28,HMF004229A,1,1
# 29,HMF004083A,0,1
# 30,HMF003821A,1,0
# 31,HMF002026A,0,1
# 32,HMF001506A,1,1
# 33,HMF001411A,1,1
# 34,HMF000685A,1,1
# 35,HMF002672A,0,1
# 36,HMF003990A,0,1

plot3 <- ggvenn(data_tn,
  fill_color = c("#0073C2FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4, show_percentage=FALSE
  ) + labs(title= "True Negatives") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
      strip.background = ggplot2::element_rect(colour = 'white', fill = 'white'),
      panel.background = element_rect(fill='white', colour='white'),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank())

data_fn <- read.csv("./Venn_predictions_fn.csv", header = TRUE, sep = ",")
data_fn$genomic = as.logical(data_fn$genomic)
data_fn$transcriptomic = as.logical(data_fn$transcriptomic)
# ,sample_ids,genomic,transcriptomic
# 0,HMF001410A,0,1
# 1,HMF001039A,1,1
# 2,HMF000574A,1,0
# 3,HMF000679A,1,0
# 4,HMF000050A,0,1
# 5,HMF001189A,0,1
# 6,HMF004001A,1,0
# 7,HMF004221A,0,1
# 8,HMF004100A,1,1
# 9,HMF004355A,0,1
# 10,HMF001378B,0,1
# 11,HMF000377A,0,1
# 12,HMF000822A,0,1
# 13,HMF000986A,0,1
# 14,HMF004057A,1,0
# 15,HMF004363A,0,1
# 16,HMF003756A,1,1
# 17,HMF004387A,0,1
# 18,HMF000759A,0,1
# 19,HMF000053A,1,1

plot4 <- ggvenn(data_fn,
  fill_color = c("#0073C2FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4, show_percentage=FALSE
  ) + labs(title= "False Negatives") +
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
      strip.background = ggplot2::element_rect(colour = 'white', fill = 'white'),
      panel.background = element_rect(fill='white', colour='white'),
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.ticks.x = element_blank(),
      axis.ticks.y = element_blank())

gl <- list(plot1, plot2, plot3, plot4)
grid.arrange(
  grobs = gl,
  layout_matrix = rbind(c(1, 2),
                        c(3, 4))
)

dev.off()