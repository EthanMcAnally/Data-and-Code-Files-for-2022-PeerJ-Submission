library(tidyverse)
library(ggplot2)
library(ggpubr)

data <- read.csv("/Users/ethanmcanally/Documents/Fall 2022/Research/Data Files/extraction_data")
extraction_quality_PS <- data[data$extraction == "ps", ]
extraction_quality_PS <- extraction_quality_PS[!duplicated(extraction_quality_PS$number),]
extraction_quality_BT <- data[data$extraction == "bt"&data$number != "", ]
extraction_quality_BT <- extraction_quality_BT[!duplicated(extraction_quality_BT$number),]
extraction_quality <- rbind(extraction_quality_BT, extraction_quality_PS)

#concentration violin plot
vc <- ggplot(extraction_quality, aes(x=extraction, y=extraction_concentration, fill=extraction)) +
  geom_violin(color = NA, trim=FALSE) +
  geom_boxplot(width=0.06) +
  ylab("ng/microliter") +
  xlab("Extraction Method") +
  scale_fill_manual(values = c("#A3CC7A","#C7B8E6")) +
  ylim(0,900) +
  theme(legend.position = "none",
        axis.ticks.x=element_blank(), plot.margin = margin(t=10,r= 10,l=10,b=10))
vc


#260/280 violin plot
vt <- ggplot(extraction_quality, aes(x=extraction, y=extraction_260.280, fill=extraction)) +
  geom_violin(color = NA, trim=FALSE) +
  geom_boxplot(width=0.07, show.legend =FALSE) +
  ylab("260/280 Absorbance") +
  xlab("Extraction Method") +
  scale_fill_manual(values = c("#A3CC7A","#C7B8E6")) +
  theme(legend.position = "none",
        axis.ticks.x=element_blank(), plot.margin = margin(t=10,r=10,l=10,b=10))
vt


#combined violin plots for figure 1
violin<- ggarrange(vc,vt, ncol=2,nrow=1) +
  theme(aspect.ratio = 0.5) 
violin

