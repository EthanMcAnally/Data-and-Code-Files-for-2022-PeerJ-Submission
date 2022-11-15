

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


###imports sample data
sample_data <- read.csv("*********/sample_prep_data_WB-FC.csv") #import csv
barcoded <- sample_data[sample_data$barcode != 00,] #remove unbarcoded samples
R_sample_data <- barcoded[barcoded$primer == "rbcl",] #makes rbcl df
S_sample_data <- barcoded[barcoded$primer == "18s",] #makes 18s df
R_sample_data <- R_sample_data[order(R_sample_data$barcode),] #orders df by barcode number
S_sample_data <- S_sample_data[order(S_sample_data$barcode),] #orders df by barcode number

###This whole block formats the raw 18s tax free csv to cleaned relative OTU abundance data from both flow cells clustered together (for figure 3)
S_data <- read.csv("********/tax_free_COMB_S.csv") #import data
S_data[,2:65] = apply(S_data[,2:65], 2, function(x) as.numeric(x)) #converts abundance columns to numeric
rownames(S_data) <- S_data[,1] #assigns first column to rownames
S_data <- S_data[,-1]#removes first column
S_data <- t(S_data)#transposes data
S_sample_names <- rownames(S_data) #make a variable of sample_names for later use
S_data <- make_relative(S_data) #convert to relative as opposed to numeric abundance
S_data <- as.data.frame(S_data) #convert back to data frame (better for manipulation)
S_data<- cbind(S_sample_names, S_data) #puts column of sample_names back in
S_data <- S_data[order(S_data$S_sample_names),]#orders data by sample names (barcode number)
S_data <- S_data[,-1] #removes column of sample names to give unfiltered relative abundance

##selects data from only the most abundant 150 taxa to match ordination processing
other <- t(S_data) #creates new matrix of transposed S_data
other <- as.data.frame(other) #convert to df
other$sums <- rowSums(other) #makes a new column of sums of overall abundance for each taxa
other$sums <- as.numeric(other$sums) #converts sums column to numeric
other <- other[order(other$sums),] #orders data by increasing taxa abundance
other <- other[122830:122979,] #only the top 150 most abundant taxa are retained to match ordination

#makes data relative again for reproducibility analysis
S_data <- t(other) #transposes data back to origional form
S_data<- as.data.frame(S_data) #converts to df
S_data <- S_data[1:64,] #removes sums
S_data<- as.matrix(S_data) #converts to matrix
S_data <- make_relative(S_data) #resets data to relative abundance
S_data <- as.data.frame(S_data) #convert back to data frame (better for manipulation)
S_data <- S_data[order(rownames(S_data)),] #orders S_data by sample names (number) barcode number
write.csv(R_data,"/Users/mcanally23/Downloads/tax_free_18S.csv", row.names = FALSE)

###This file was then modified to the "Diatom_taxfree_18S.csv" file for reproducibility analysis


#This whole block formats the raw rbcL tax free csv to cleaned relative OTU abundance data from the first flow cell (for figure 4)
R_data <- read.csv("******/tax_free_WB_R.csv") #import data
R_data[,2:31] = apply(R_data[,2:31], 2, function(x) as.numeric(x))
rownames(R_data) <- R_data[,1]
R_data <- R_data[,-1]
R_data <- t(R_data)
R_sample_names <- rownames(R_data) #make a variable of sample_names for later use
R_data <- make_relative(R_data) #convert to relative as opposed to numeric abundance
R_data <- as.data.frame(R_data) #convert back to data frame (better for manipulation)
R_data<- cbind(R_sample_names, R_data) #puts column of sample_names back in
R_data <- R_data[order(R_data$R_sample_names),]#orders data by barcode number
R_data <- R_data[,-1]
com_R <- R_data

other <- t(R_data)
other <- as.data.frame(other)
other$sums <- rowSums(other)
other$sums <- as.numeric(other$sums)
other <- other[order(other$sums),]
other <- other[70967:71116,]
R_data <- t(other)
R_data<- as.data.frame(R_data)
R_data <- R_data[1:30,]
R_data<- as.matrix(R_data)
R_data <- make_relative(R_data) #resets data to relative abundance (otherwise would not add to 100%)
R_data <- as.data.frame(R_data) #convert back to data frame (better for manipulation)
R_data <- R_data[order(rownames(R_data)),]
write.csv(R_data,"/Users/mcanally23/Downloads/tax_free_rbcL.csv", row.names = FALSE)

###This file was then modified to the "Diatom_taxfree_rbcL.csv" file for NMDS analysis


