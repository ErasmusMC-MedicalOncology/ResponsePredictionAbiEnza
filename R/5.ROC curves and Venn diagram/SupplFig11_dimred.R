library(ggplot2)
library(RColorBrewer)
library(patchwork)
require(gridExtra)
library(grid)
library(cowplot)
library(stringr)

#### TP (true positives), FP (false positives), TN (true negatives) and FN (false negatives) can be extracted from
#### the Suppl. Table 1 / Source Data to calculate AUC and ROC curves
#### predictions for models: "Predictions - HMF" sheet
#### ground truth labels: "Sample information" sheet ("Responder")

cols_plot3 <- c("10 ICs AUC=0.53" = "#7ac167", "15 ICs AUC=0.53" = "#3f8f38", "20 ICs AUC=0.74" = "#B3DE69",
"25 ICs AUC=0.74" = "#778a35", "30 ICs AUC=0.70" = "#66C2A5", "35 ICs AUC=0.65" = "#1B9E77", "40 ICs AUC=0.76" = "#76b947", "50 ICs AUC=0.74" = "#a0de7c")
# blue
cols_plot1 <- c("10 sparse PCs AUC=0.38" = "#8face7", "15 sparse PCs AUC=0.38" = "#6087d5", "20 sparse PCs AUC=0.39" = "#3965bd",
"25 sparse PCs AUC=0.42" = "#1b3cdf", "30 sparse PCs AUC=0.40" = "#1b8cdf", "35 sparse PCs AUC=0.41" = "#448abd", "40 sparse PCs AUC=0.51" = "#4dadf3", "50 sparse PCs AUC=0.58" = "#5fbdde")
# red
cols_plot2 <- c("10 PCs AUC=0.33" = "#e15804", "15 PCs AUC=0.33" = "#ad5926", "20 PCs AUC=0.63" = "#e83809", "25 PCs AUC=0.56" = "#e16949",
"30 PCs AUC=0.54" = "#be401e",
          "35 PCs AUC=0.63" = "#f52d19", "40 PCs AUC=0.65" = "#f76556", "50 PCs AUC=0.66" = "#f79156")


pc_data <- read.csv("./ROC_R_input_Suppl_Fig11_dimred_PCA.csv", header = TRUE, sep = ",")
ic_data <- read.csv("./ROC_R_input_Suppl_Fig11_dimred_ICA.csv", header = TRUE, sep = ",")
spc_data <- read.csv("./ROC_R_input_Suppl_Fig11_dimred_sPCA.csv", header = TRUE, sep = ",")
spc_data <- subset(pc_data, select = -c(X))


plot3 <- ggplot(ic_data,aes(fpr,tpr,color=model))+geom_line(size = 2, alpha = 0.7) +
      labs(fill='independent components', title = "c", x = " ", y=" ") + # + theme_bw()
      theme(
      # t, r, b, l
#       plot.margin = unit(c(1,0,1,1), "cm"),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position="none",
      text = ggplot2::element_text(size=12, family='Helvetica', face = 'bold'),
      axis.text.y = element_text(size=10, family='Helvetica', face = 'bold', margin = margin(t = 0, r = 10, b = 0, l = 10)),
      axis.text.x = element_text(size=10, family='Helvetica', face = 'bold', margin = margin(t = 10, r = 0, b = 10, l = 0)),
      strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
      panel.background = element_rect(fill='white', colour='black'),
      panel.grid.major = ggplot2::element_line(colour = 'grey', linetype = 'dotted'),
      panel.grid.minor = ggplot2::element_blank()) + scale_colour_manual(values = cols_plot3, aesthetics = c("colour", "fill"))

plot2 <- ggplot(pc_data,aes(fpr,tpr,color=model))+geom_line(size = 2, alpha = 0.7) +
      labs(fill='principal components', title = "b", x = "False Positive Rate", y=" ") +
      theme(
      # t, r, b, l
#       plot.margin = unit(c(1,0,1,0), "cm"),
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.position="none",
      text = ggplot2::element_text(size=12, family='Helvetica', face = 'bold'),
      axis.text.x = element_text(size=10, family='Helvetica', face = 'bold', margin = margin(t = 10, r = 10, b = 10, l = 0)),
      axis.text.y = element_text(size=10, family='Helvetica', face = 'bold', margin = margin(t = 0, r = 10, b = 0, l = 10)),
      strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
      panel.background = element_rect(fill='white', colour='black'),
      panel.grid.major = ggplot2::element_line(colour = 'grey', linetype = 'dotted'),
      panel.grid.minor = ggplot2::element_blank()) + scale_colour_manual(values = cols_plot2, aesthetics = c("colour", "fill"))

plot1 <- ggplot(spc_data,aes(fpr,tpr,color=model))+geom_line(size = 2, alpha = 0.7) +
      labs(fill='sparse principal components', title = "a", x = " ", y = "True Positive Rate") + # + theme_bw()
      theme(
      legend.position="none",
      text = ggplot2::element_text(size=12, family='Helvetica', face = 'bold'),
      axis.text.x = element_text(size=10, family='Helvetica', face = 'bold', margin = margin(t = 10, r = 0, b = 10, l = 0)),
      axis.text.y = element_text(size=10, family='Helvetica', face = 'bold', margin = margin(t = 0, r = 10, b = 0, l = 10)),
      strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
      panel.background = element_rect(fill='white', colour='black'),
      panel.grid.major = ggplot2::element_line(colour = 'grey', linetype = 'dotted'),
      panel.grid.minor = ggplot2::element_blank()) + guides(fill=guide_legend(ncol=2)) + guides(color=guide_legend(ncol=2)) + scale_colour_manual(values = cols_plot1)


png(file="./Suppl_Fig_11_dimred.png", width=20, height = 9, units = "in", res= 800, pointsize = 4)

gl <- list(plot1, plot2, plot3)
# t, r, b, l
margin = theme(plot.margin = unit(c(1,1.5,1,0), "cm"))
# grid.arrange(grobs = lapply(gl, "+", margin))
grid.arrange(
  grobs = lapply(gl, "+", margin),
  widths = c(1.45, 1.40, 1.40),
  heights = c(1, 1, 1),
  layout_matrix = rbind(c(1, 2, 3),
                        c(1, 2, 3),
                        c(1, 2, 3))
)
dev.off()
