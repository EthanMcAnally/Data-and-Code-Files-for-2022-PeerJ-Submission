
library(tidyverse)
library(ggplot2)
library(funrar)
library(dplyr)
library(readr)
library(reshape2)
library(vegan)
library(ggthemes)
library(extrafont)
library(scales)
library(ggpubr)
library(gridExtra)
library(rstatix) 
library(emmeans) 
library(MASS)
library(pals)
library(Polychrome)

#imports non-genetic data on samples
sample_data <- read.csv("/Users/ethanmcanally/Documents/Fall 2022/Research/Data Files/sample_prep_data_WB-FC.csv") #import csv
barcoded <- sample_data[sample_data$barcode != 00,] #remove unbarcoded samples
R_sample_data <- barcoded[barcoded$primer == "rbcl",] #makes rbcl df
S_sample_data <- barcoded[barcoded$primer == "18s",] #makes 18s df
R_sample_data <- R_sample_data[order(R_sample_data$barcode),] #orders df by barcode number
S_sample_data <- S_sample_data[order(S_sample_data$barcode),] #orders df by barcode number

##imports and formats 18S relative tax abundance data from Qiime2
#WB=code for first flow cell run S=18S
WB_S_data <- read.csv("/Users/ethanmcanally/Documents/Fall 2022/Research/Data Files/tax_WB_S.csv") #import data
WB_S_sample_names <- WB_S_data[,1] #make a variable of sample_names for later use
WB_S_data <- WB_S_data[,-1] #remove column of sample names, so all values are numeric
wildcard <- "*Bacillariophyta*"
vector <- colnames(WB_S_data)
sum(WB_S_data)
keepers <- grep(wildcard, vector, value = TRUE)
WB_S_data <- WB_S_data[c(keepers)]
sum(WB_S_data)
vector <- character()
for (x in colnames(WB_S_data)) {
  x <- gsub(".", ";", x, fixed = TRUE)
  x <- gsub(";__", "", x, fixed = TRUE)
  x <- gsub(".*;", "", x,)
  vector <- c(vector,x)
}
colnames(WB_S_data) <- vector
WB_S_data <- data.matrix(WB_S_data) #convert to matrix (necessary for next step)
WB_S_data <- make_relative(WB_S_data) #convert to relative as opposed to numeric abundance
other <- WB_S_data
other[other > 0.01] <- 0
other <- rowSums(other)
WB_S_data[WB_S_data < 0.01] <- 0 #removes data on species composing less than 1.0% of a sample
WB_S_data <- cbind(WB_S_data,other)
WB_S_data <- make_relative(WB_S_data) #resets data to relative abundance (otherwise would not add to 100%)
WB_S_data <- as.data.frame(WB_S_data) #convert back to data frame (better for manipulation)
WB_S_data <- WB_S_data[, colSums(WB_S_data != 0) > 0]
WB_S_data<- cbind(WB_S_sample_names, WB_S_data) #puts column of sample_names back in
WB_S_data <- WB_S_data[order(WB_S_data$WB_S_sample_names),]#orders data by barcode number
WB_S_data <- cbind(S_sample_data, WB_S_data)#adds columns with non-genetic data
WB_S_data<- WB_S_data[,-12] #removes old sample names column
WB_S_data <- melt(WB_S_data, id.vars = 1:11, variable.name = "Phyla") #converts to long form for ggplotting
WB_S_data <- WB_S_data[WB_S_data$value != 0,] #removes taxa that are no longer in any sample (after 1% cutoff)
WB_S_data$barcode <- as.factor(WB_S_data$barcode)#assigns class of barcode column to factor for better figure
WB_S_data <- WB_S_data[order(WB_S_data$value),]#orders data by abundance for better looking figure

##imports and formats rbcL relative tax abundance data from Qiime2
#WB=code for first flow cell run R=rbcL
WB_R_data <- read.csv("/Users/ethanmcanally/Documents/Fall 2022/Research/Data Files/tax_WB_R.csv") #import data
WB_R_sample_names <- WB_R_data[,1] #make a variable of sample_names for later use
WB_R_data <- WB_R_data[,-1] #remove column of sample names, so all values are numeric
wildcard <- "*Bacillariophyta*"
vector <- colnames(WB_R_data)
sum(WB_R_data)
keepers <- grep(wildcard, vector, value = TRUE)
WB_R_data <- WB_R_data[c(keepers)]
sum(WB_R_data)
vector <- character()
for (x in colnames(WB_R_data)) {
  x <- gsub(".", ";", x, fixed = TRUE)
  x <- gsub(";__", "", x, fixed = TRUE)
  x <- gsub(".*;", "", x,)
  vector <- c(vector,x)
}
colnames(WB_R_data) <- vector
WB_R_data <- data.matrix(WB_R_data) #convert to matrix (necessary for next step)
WB_R_data <- make_relative(WB_R_data) #convert to relative as opposed to numeric abundance
other <- WB_R_data
other[other > 0.01] <- 0
other <- rowSums(other)
WB_R_data[WB_R_data < 0.01] <- 0 #removes data on species composing less than 1.0% of a sample
WB_R_data <- cbind(WB_R_data,other)
WB_R_data <- make_relative(WB_R_data) #resets data to relative abundance (otherwise would not add to 100%)
WB_R_data <- as.data.frame(WB_R_data) #convert back to data frame (better for manipulation)
WB_R_data <- WB_R_data[, colSums(WB_R_data != 0) > 0]
WB_R_data<- cbind(WB_R_sample_names, WB_R_data) #puts column of sample_names back in
WB_R_data <- WB_R_data[order(WB_R_data$WB_R_sample_names),]#orders data by barcode number
WB_R_data <- cbind(R_sample_data, WB_R_data)#adds columns with non-genetic data
WB_R_data<- WB_R_data[,-12] #removes old sample names column
WB_R_data <- melt(WB_R_data, id.vars = 1:11, variable.name = "Phyla") #converts to long form for ggplotting
WB_R_data <- WB_R_data[WB_R_data$value != 0,] #removes taxa that are no longer in any sample (after 1% cutoff)
WB_R_data$barcode <- as.factor(WB_R_data$barcode)#assigns class of barcode column to factor for better figure
WB_R_data <- WB_R_data[order(WB_R_data$value),]#orders data by abundance for better looking figure

###imports microscope based data
M_data <- read.csv("/Users/ethanmcanally/Documents/Fall 2022/Research/Data Files/Microscopic_genus_abundance.csv") #import csv
rownames(M_data) <- M_data[,1] #assigns sample names column to row names
M_data <- M_data[,-1] #remove column of sample names, so all values are numeric
M_data[M_data < 0.01] <- 0 #removes data on species composing less than 1.0% of a sample
M_data <- data.matrix(M_data) #convert to matrix (necessary for next step)
M_data <- make_relative(M_data) #resets data to relative abundance (otherwise would not add to 100%)
M_data <- as.data.frame(M_data) #convert back to data frame (better for manipulation)
M_data <- M_data[order(rownames(M_data)),]#orders data by barcode to match non-genetic data
M_data$sample <- rownames(M_data)#adds column with barcode names
#takes subset of data as only the sample relevant for ecological analysis
M_data <- M_data[M_data$sample=="1"|M_data$sample=="2"|M_data$sample=="6"|M_data$sample=="7"|M_data$sample=="11"|M_data$sample=="12"|M_data$sample=="16"|M_data$sample=="17",]
M_data <- M_data[order(M_data$sample),]#orders data by barcode to match non-genetic data
M_sample_data <- S_sample_data[S_sample_data$barcode=="40"|S_sample_data$barcode=="28"|S_sample_data$barcode=="77"|S_sample_data$barcode=="65"|S_sample_data$barcode=="17"|S_sample_data$barcode=="5"|S_sample_data$barcode=="54"|S_sample_data$barcode=="42",]
M_sample_data <- M_sample_data[order(M_sample_data$number),]#orders data by barcode to match non-genetic data
M_data <- M_data[,-24]
M_data <- cbind(M_sample_data,M_data)
M_data <- melt(M_data, id.vars = 1:11, variable.name = "Taxon")
M_data <- M_data[M_data$value != 0,]
M_list <- as.list(rep("Microscopic", 68))
M_data$primer <- as.character(M_list)
M_data$value <- as.numeric(M_data$value)
M_data$number <- as.character(M_data$number)

#######prepares data for primer comparison figure (figure 2b)
primers <- rbind(WB_S_data,WB_R_data)
names(primers)[names(primers) == "Phyla"] <- "Taxon"
primers <- primers[primers$barcode=="40"|primers$barcode=="28"|primers$barcode=="77"|primers$barcode=="65"|primers$barcode=="17"|primers$barcode=="5"|primers$barcode=="54"|primers$barcode=="42"|primers$barcode=="86"|primers$barcode=="74"|primers$barcode=="26"|primers$barcode=="14"|primers$barcode=="63"|primers$barcode=="51"|primers$barcode=="3"|primers$barcode=="88",]
primers <- rbind(primers,M_data)
primers <- primers[primers$value != 0,]
primers$Taxon <- as.character(primers$Taxon)
primers$sort <- primers$Taxon
#changes letter codes to order by taxonomy, will correct in pdf
#I realize this is inelegant :)
primers$sort[primers$sort == "f__Achnanthidiaceae"] <- "h__Achnanthidiaceae"
primers$sort[primers$sort == "f__Bacillariaceae"] <- "h__Bacillariaceae"
primers$sort[primers$sort == "f__Fragilariaceae"] <- "h__Fragilariaceae"
primers$sort[primers$sort == "f__Gomphonemataceae"] <- "h__Gomphonemataceae"
primers$sort[primers$sort == "other"] <- "q__other"
primers <- primers[order(primers$sort),]
primers$number <- as.numeric(primers$number)
primers <- primers[order(primers$number),]
primers$number <- as.factor(primers$number)
primers[primers == "18s"] <- "amplicon_18s"
primers[primers == "rbcl"] <- "amplicon_rbcl"
primer_labs <- c("18s Genetic Data", "rbcL Genetic Data")
names(primer_labs) <- c("18s Genetic Data", "rbcL Genetic Data")
palette_p <- c("#990F26","#B33E52","#CC7A88","#E6B8BF","#99600F","#B3823E","#CCAA7A","#E6D2B8","#54990F","#78B33E","#A3CC7A","#CFE6B8","#0F8299","#3E9FB3","#7ABECC","#B8DEE6","#3D0F99","#653EB3","#967ACC","#C7B8E6","#333333","#666666","#999999","#CBCBCB", "#DEDEDE", "#FFFFFF" )

tax_amplification_figure<-ggplot(primers, aes(fill=sort, y=value, x=number)) + 
  geom_bar(position="stack", stat="identity") +
  facet_grid(cols = vars(primer)) +
  xlab("Sample") +
  ylab("Relative Abundance") +
  scale_fill_manual(values = palette_p) +
  guides(fill=guide_legend(ncol =1)) +
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5, size = 20), plot.margin = margin(t=0,r=0,l=10,b=0), legend.text = element_text(size = 8), legend.key.size = unit(0.9,"line"), legend.key.width = unit(0.5, 'cm'))
tax_amplification_figure

#####imports and processes data for relative Amplfication Success Figure (panel 2a)
#import data file
data <- read.csv("/Users/ethanmcanally/Documents/Fall 2022/Research/Data Files/sample_prep_data _WB-FC.csv")
#all 18s, ps samples that were successfully amplified
SSAP <- data[data$lost!="pcr"&data$lost!="pcr_c"&data$lost!="extraction"&data$lost!="extraction_left" & data$extraction == "ps" & data$primer == "18s",]
#all 18s, ps samples that failed amplification
SSAF <- data[data$lost!="amplicons_lost"&data$lost!="barcoded"&data$lost!="misbarcoded"&data$lost!="extraction"&data$lost!="extraction_left" & data$extraction == "ps" & data$primer == "18s",]
#repeating this for pass and fail of all four combinations
SRAP <- data[data$lost!="pcr"&data$lost!="pcr_c"&data$lost!="extraction"&data$lost!="extraction_left" & data$extraction == "ps" & data$primer == "rbcl",]
SRAF <- data[data$lost!="amplicons_lost"&data$lost!="barcoded"&data$lost!="misbarcoded"&data$lost!="extraction"&data$lost!="extraction_left" & data$extraction == "ps" & data$primer == "rbcl",]
BSAP <- data[data$lost!="pcr"&data$lost!="pcr_c"&data$lost!="extraction"&data$lost!="extraction_left" & data$extraction == "bt" & data$primer == "18s",]
BSAF <- data[data$lost!="amplicons_lost"&data$lost!="barcoded"&data$lost!="misbarcoded"&data$lost!="extraction"&data$lost!="extraction_left" & data$extraction == "bt" & data$primer == "18s",]
BRAP <- data[data$lost!="pcr"&data$lost!="pcr_c"&data$lost!="extraction"&data$lost!="extraction_left"& data$extraction == "bt" & data$primer == "rbcl",]
BRAF <- data[data$lost!="amplicons_lost"&data$lost!="barcoded"&data$lost!="misbarcoded"&data$lost!="extraction"&data$lost!="extraction_left" & data$extraction == "bt" & data$primer == "rbcl",]

#creating df of success/failure outcomes
extraction <- c("PS","PS","PS","PS","BT","BT","BT","BT")
primer <- c("18S","18S","rbcL","rbcL","18S","18S","rbcL","rbcL")
outcome <- c("success","failure","success","failure","success","failure","success","failure")
sample_count <- c(nrow(SSAP),nrow(SSAF),nrow(SRAP),nrow(SRAF),nrow(BSAP),nrow(BSAF),nrow(BRAP),nrow(BRAF))

amplification_df <- data.frame(extraction, primer, outcome, sample_count)

print (amplification_df)

#adds row of calculated relative sucesses and failures (must be a better way, but this was quick)
relative <- c(1,0,1,0,0.667,0.333,0.333,0.667)
amplification_df <- cbind(amplification_df,relative)

Palette <- c("#990F26","#78B33E")

amplification_figure<- ggplot(data=amplification_df, aes(x=extraction, y=relative, fill=outcome)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values = Palette) +
  scale_x_discrete(name = "Extraction Method") +
  scale_y_continuous(name = "Frequency of Amplification Success") +
  facet_grid(cols = vars(primer)) +
  theme(text = element_text(size = 10), legend.position = "none" , strip.background = element_blank(), plot.margin = margin(t=0,r=0,l=0,b=0))
amplification_figure


combined_amplification_figure <- ggarrange(amplification_figure,tax_amplification_figure, ncol=2,nrow=1,widths = c(.3, 1)) +
  theme(plot.margin = margin(10,0,10,10)) 
combined_amplification_figure


#environmental comparison figure
env <- R_A_data
env <- env[env$barcode=="2"|env$barcode=="14"|env$barcode=="9"|env$barcode=="8"|env$barcode=="20"|env$barcode=="32"|env$barcode=="74"|env$barcode=="86"|env$barcode=="3"|env$barcode=="15"|env$barcode=="27"|env$barcode=="39"|env$barcode=="51"|env$barcode=="63"|env$barcode=="75",]
write.csv(env,"/Users/mcanally23/Downloads/environmental_relative_taxa_abundance.csv", row.names = FALSE)
####this data was modified for NDMS analysis to make the "Diatom_taxbased_18S.csv" and "Diatom_taxbased_rbcL.csv" files uploaded in the repository

