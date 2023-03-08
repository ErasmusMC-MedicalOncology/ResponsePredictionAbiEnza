library(ggplot2)
library(RColorBrewer)
library(patchwork)
require(gridExtra)
library(grid)
library(cowplot)
library(stringr)

cols_plot5 <- c("genomics AUC=0.76" = "#E41A1C", "genomics + ARSI AUC=0.81" = "#297e9f", "transcriptomics AUC=0.76" = "#778a35", "transcriptomics + ARSI AUC=0.82" = "#0fba0c", "transcriptomics + genomics (avg. ensemble) AUC=0.81" = "#984EA3",
"transcriptomics + genomics + ARSI AUC=0.84" = "#b7e78f")
cols_plot4 <- c("genomics" = "#E41A1C", "genomics + ARSI" = "#297e9f", "transcriptomics" = "#3f8f38", "transcriptomics + ARSI" = "#0fba0c", "avg. ensemble" = "#984EA3", "transcriptomics + genomics + ARSI" = "#b7e78f")


cols_plot1 <- c("10 ICs AUC=0.53" = "#7ac167", "15 ICs AUC=0.53" = "#3f8f38", "20 ICs AUC=0.74" = "#B3DE69",
"25 ICs AUC=0.74" = "#778a35", "30 ICs AUC=0.70" = "#66C2A5", "35 ICs AUC=0.65" = "#1B9E77", "40 ICs AUC=0.76" = "#76b947", "50 ICs AUC=0.74" = "#a0de7c", "genomics AUC=0.76" = "#E41A1C")
cols_plot2 <- c("multi-model avg. ensemble AUC=0.75" = "#BC80BD", "avg. ensemble AUC=0.81" = "#984EA3", "stacking classifier AUC=0.75" = "#e23dad", "bagging classifier AUC=0.75" = "#93518d")
cols_plot3 <- c("genomics + ARSI + chemo AUC=0.80" = "#3b99b8", "genomics + ARSI AUC=0.81" = "#297e9f", "genomics + chemo AUC=0.76" = "#4fb5d0", "genomics + n.t.l. AUC=0.77" = "#64d1e8", "genomics AUC=0.76" = "#E41A1C", "ARSI + chemo AUC=0.61" = "#808080", "transcriptomics + ARSI AUC=0.82" = "#0fba0c", "transcriptomics + genomics + ARSI AUC=0.84" = "#b7e78f")

sdata <- read.csv("./ROC_R_best_models.csv", header = TRUE, sep = ",")
sdata <- subset(sdata, select = -c(X))



sdata$model <- factor(sdata$model, levels = c("genomics AUC=0.76", "genomics + ARSI AUC=0.81",
"transcriptomics AUC=0.76", "transcriptomics + ARSI AUC=0.82", "transcriptomics + genomics (avg. ensemble) AUC=0.81",
"transcriptomics + genomics + ARSI AUC=0.84"))

plot5 <- ggplot(sdata,aes(fpr,tpr,color=model))+geom_line(size = 2, alpha = 0.7) +
      labs(title= "e",
           x = "False Positive Rate",
           y = "True Positive Rate") + # + theme_bw()
      theme(
      legend.background = element_blank(),
      legend.box.background = element_rect(fill = alpha('white', 0.2)),
      legend.position = c(.70, .2),
      legend.direction = 'vertical',

      legend.text = element_text(size=8, family='Helvetica', face = 'bold'),
      legend.key=element_blank(),
      text = ggplot2::element_text(size=12, family='Helvetica', face = 'bold'),
      axis.text.x = element_text(size=12, family='Helvetica', face = 'bold', margin = margin(t = 10, r = 0, b = 10, l = 0)),
      axis.text.y = element_text(size=12, family='Helvetica', face = 'bold', margin = margin(t = 0, r = 10, b = 0, l = 10)),
      strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
      panel.background = element_rect(fill='white', colour='black'),
      panel.grid.major = ggplot2::element_line(colour = 'grey', linetype = 'dotted'),
      panel.grid.minor = ggplot2::element_blank()) + scale_colour_manual(values = cols_plot5)+ guides(color = guide_legend(override.aes = list(size = 5)))

data <- read.csv("/shuffled_LOOCV_experiments.csv", header = TRUE, sep = ",")
data$experiment <- factor(data$experiment, levels = c("genomics", "genomics + ARSI", "transcriptomics", "transcriptomics + ARSI", "avg. ensemble", "transcriptomics + genomics + ARSI"))
plot4 <- ggplot(data, aes(x=experiment, y=value, fill=experiment)) + geom_boxplot(width=0.5,lwd=1) +
      labs(title = "d", y = 'AUC in shuffled label experiment', x = '') +
      theme(
      legend.background = element_blank(),
      legend.box.background = element_rect(fill = alpha('white', 0.2)),
      legend.position = "none",
      legend.direction = 'vertical',
      legend.text = element_text(size=8, family='Helvetica', face = 'bold'),
      legend.key.size = unit(0.8, 'cm'),
      text = ggplot2::element_text(size=12, family='Helvetica', face = 'bold'),
      axis.text.x = element_text(size=12, family='Helvetica', face = 'bold', margin = margin(t = 10, r = 0, b = 10, l = 0)),
      axis.text.y = element_text(size=12, family='Helvetica', face = 'bold', margin = margin(t = 0, r = 10, b = 0, l = 10)),
      strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
      panel.background = element_rect(fill='white', colour='black'),
      panel.grid.major = ggplot2::element_line(colour = 'grey', linetype = 'dotted'),
      panel.grid.minor = ggplot2::element_blank()) +  scale_colour_manual(values = cols_plot4, aesthetics = c("colour", "fill"))
plot4 + scale_x_discrete(labels = function(x) str_wrap(x, width = 5))

ensembledata <- read.csv("./ROC_R_input_ensembles.csv", header = TRUE, sep = ",")
cdata <- read.csv("./ROC_all_ONLY_ccov_standalone_and_with_genomics.csv", header = TRUE, sep = ",")
pc_data <- read.csv("./ROC_R_input_Fig5_a.csv", header = TRUE, sep = ",")
pc_data <- subset(pc_data, select = -c(X))


plot3 <- ggplot(cdata,aes(fpr,tpr,color=model))+geom_line(size = 1, alpha = 0.7) +
      labs(fill='clinical covariate', title = "c", x = " ", y=" ") + # + theme_bw()
      theme(
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.background = element_blank(),
      legend.position = c(0.55, 0.35),
      legend.direction = 'vertical',
      legend.title = element_blank(),
      legend.text = element_text(size=8, family='Helvetica', face = 'bold'),
      legend.key=element_blank(),
      text = ggplot2::element_text(size=12, family='Helvetica', face = 'bold'),
      axis.text.y = element_text(size=10, family='Helvetica', face = 'bold', margin = margin(t = 0, r = 10, b = 0, l = 10)),
      axis.text.x = element_text(size=10, family='Helvetica', face = 'bold', margin = margin(t = 10, r = 0, b = 10, l = 0)),
      strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
      panel.background = element_rect(fill='white', colour='black'),
      panel.grid.major = ggplot2::element_line(colour = 'grey', linetype = 'dotted'),
      panel.grid.minor = ggplot2::element_blank()) + scale_colour_manual(values = cols_plot3, aesthetics = c("colour", "fill"))

plot2 <- ggplot(ensembledata,aes(fpr,tpr,color=model))+geom_line(size = 1, alpha = 0.7) +
      labs(fill='ensemble models', title = "b", x = "False Positive Rate", y=" ") +
      theme(
      axis.ticks.y = element_blank(),
      axis.title.y = element_blank(),
      legend.background = element_blank(),
      legend.position = c(0.55, 0.6),
      legend.direction = 'vertical',
      legend.title = element_blank(),
      legend.text = element_text(size=8, family='Helvetica', face = 'bold'),
      legend.key=element_blank(),
      text = ggplot2::element_text(size=12, family='Helvetica', face = 'bold'),
      axis.text.x = element_text(size=10, family='Helvetica', face = 'bold', margin = margin(t = 10, r = 10, b = 10, l = 0)),
      axis.text.y = element_text(size=10, family='Helvetica', face = 'bold', margin = margin(t = 0, r = 10, b = 0, l = 10)),
      strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
      panel.background = element_rect(fill='white', colour='black'),
      panel.grid.major = ggplot2::element_line(colour = 'grey', linetype = 'dotted'),
      panel.grid.minor = ggplot2::element_blank()) + scale_colour_manual(values = cols_plot2, aesthetics = c("colour", "fill"))

plot1 <- ggplot(pc_data,aes(fpr,tpr,color=model))+geom_line(size = 1, alpha = 0.7) +
      labs(fill='sparse PCs', title = "a", x = " ", y = "True Positive Rate") + # + theme_bw()
      theme(
      legend.background = element_blank(),
      legend.position = c(0.6, 0.35),
      legend.direction = 'vertical',
      legend.title = element_blank(),
      legend.text = element_text(size=8, family='Helvetica', face = 'bold'),
      legend.key=element_blank(),
      text = ggplot2::element_text(size=12, family='Helvetica', face = 'bold'),
      axis.text.x = element_text(size=10, family='Helvetica', face = 'bold', margin = margin(t = 10, r = 0, b = 10, l = 0)),
      axis.text.y = element_text(size=10, family='Helvetica', face = 'bold', margin = margin(t = 0, r = 10, b = 0, l = 10)),
      strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
      panel.background = element_rect(fill='white', colour='black'),
      panel.grid.major = ggplot2::element_line(colour = 'grey', linetype = 'dotted'),
      panel.grid.minor = ggplot2::element_blank()) + guides(fill=guide_legend(ncol=2)) + guides(color=guide_legend(ncol=2)) + scale_colour_manual(values = cols_plot1)


legend_p1 <- get_legend(plot1)
legend_p2 <- get_legend(plot2)
legend_p3 <- get_legend(plot3)
legend_p1$grobs[[1]]$widths[4] <- unit(3.0, "cm")

# remove legend
bp1_wl <- plot1 + theme(legend.position='none')
bp2_wl <- plot2 + theme(legend.position='none')
bp3_wl <- plot3 + theme(legend.position='none')

plot1_3 <- plot_grid(
  bp1_wl, bp2_wl, bp3_wl, legend_p1, legend_p2, legend_p3,
  ncol = 3, nrow = 2, rel_heights = c(6, 2), rel_widths = c(6.3, 6, 6)
)

png(file="./Fig_5_panel.png", width=20, height = 9, units = "in", res= 1200, pointsize = 4)

gl <- list(plot1_3, plot4, plot5)
margin = theme(plot.margin = unit(c(1,1.5,1,0), "cm"))
grid.arrange(
  grobs = lapply(gl, "+", margin),
  widths = c(1.45, 1.40, 1.40, 3),
  heights = c(1, 0.40, 1),
  layout_matrix = rbind(c(1, 1, 1, 3),
                        c(1, 1, 1, 3),
                        c(2, 2, 2, 3))
)
dev.off()
