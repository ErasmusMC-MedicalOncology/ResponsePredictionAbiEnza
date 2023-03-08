library(ggplot2)
library(RColorBrewer)
library(patchwork)
require(gridExtra)
library(grid)
library(cowplot)

cols_plot5 <- c("all genes AUC=0.74" = "#297e9f", "filtered AUC=0.74" = "#4fb5d0")

sdata <- read.csv("./ROC_R_input_filt_vs_unfilt.csv", header = TRUE, sep = ",")
sdata <- subset(sdata, select = -c(X))


sdata$model <- factor(sdata$model, levels = c("all genes AUC=0.76", "filtered AUC=0.76"))


plot2 <- ggplot(sdata,aes(fpr,tpr,color=model))+geom_line(size = 2, alpha = 0.7) +
      labs(title= "filtered/unfiltered transcriptomics (ICs 40) comparison",
           x = "False Positive Rate",
           y = "True Positive Rate") + # + theme_bw()
      theme(
      axis.ticks.y = element_blank(),
      legend.background = element_blank(),
      legend.box.background = element_rect(fill = alpha('white', 0.2)),
      legend.position = c(.52, .15),
      legend.direction = 'vertical',
      legend.title = element_blank(),
      legend.text = element_text(size=10, family='Helvetica', face = 'bold'),
      legend.key.size = unit(0.6, 'cm'),
      text = ggplot2::element_text(size=10, family='Helvetica', face = 'bold'),
      axis.text.x = element_text(size=10, family='Helvetica', face = 'bold', margin = margin(t = 10, r = 10, b = 10, l = 0)),
      axis.text.y = element_text(size=10, family='Helvetica', face = 'bold', margin = margin(t = 0, r = 10, b = 0, l = 10)),
      strip.background = ggplot2::element_rect(colour = 'black', fill = 'white'),
      panel.background = element_rect(fill='white', colour='black'),
      panel.grid.major = ggplot2::element_line(colour = 'grey', linetype = 'dotted'),
      panel.grid.minor = ggplot2::element_blank())
ggsave("./SupplFig12_filt_vs_unfilt.png")