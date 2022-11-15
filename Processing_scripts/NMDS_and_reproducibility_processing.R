#############NDS Analysis
library(ggplot2)
library(ggthemes)
library(extrafont)
library(plyr)
library(scales)
library(vegan)
library(MASS)
library(ggpubr)
library(ggplot2)
library(gridExtra) #grid arrange
library(rstatix) #pipe friendly R functions for statistics
library(emmeans) #multiple comparisons
library(gridExtra) #grid arrange
library(RColorBrewer)


##############NMDS - 2022
#MICROSCOPE
df<-read.csv("Diatom_microscope_20221030.csv")
df
colnames(df)
en<-read.csv("C:\\Users\\weilhoef\\Desktop\\eDNA\\Env_20221027.csv")
en
colnames(en)

df$HABITAT <- as.character(df$HABITAT)
df$HABITAT <- factor(df$HABITAT, levels=unique(df$HABITAT))
df$RIVER <- as.character(df$RIVER)
df$RIVER <- factor(df$RIVER, levels=unique(df$RIVER))
com = df[,5:43]
env = en[,4:12]
#convert com to a matrix
m_com = as.matrix(com)
bcd = vegdist(m_com, "bray")
bcd

#nmds code
set.seed(765)
nmds = metaMDS(m_com, distance = "bray")
nmds
en = envfit(nmds, env, permutations = 999, na.rm = TRUE)
en
plot(nmds)
plot(en)

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores = as.data.frame(scores(nmds)$sites)
data.scores$RIVER = df$RIVER #adds in column with River 
data.scores$HABITAT = df$HABITAT #adds in column with Treatment
head(data.scores)
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)

#plot
gg = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(shape = HABITAT, colour = RIVER), size = 5, alpha = 0.85) +  
  scale_shape_manual(values = c(16, 17)) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size = 0.75, alpha = 0.5, colour = "grey30",
               arrow = arrow(length=unit(0.25,"cm"), type = "closed")) +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            nudge_x = 0.085, nudge_y = -0.01,size=3, label = row.names(en_coord_cont)) + 
  theme(axis.title = element_text(size = 10, colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30")) +
  theme(legend.key = element_rect(colour = "transparent", fill = "white"), legend.position =  "none",
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
        legend.direction = "vertical", legend.box = "horizontal",
        legend.title = element_text(size=10),
        legend.text=element_text(size=8)) +
  labs(colour = "River", shape = "Habitat")
gg1 <- gg +  scale_colour_manual(values = brewer.pal(4, "Paired")[c(2,4)])
gg1

#rbcL - taxonomy based
df<-read.csv("Diatom_taxbased_rbcl.csv")
df
colnames(df)
en<-read.csv("Env_20221027.csv")
en
colnames(en)

df$HABITAT <- as.character(df$HABITAT)
df$HABITAT <- factor(df$HABITAT, levels=unique(df$HABITAT))
df$RIVER <- as.character(df$RIVER)
df$RIVER <- factor(df$RIVER, levels=unique(df$RIVER))
com = df[,5:23]
env = en[,4:12]
#convert com to a matrix
m_com = as.matrix(com)
bcd = vegdist(m_com, "bray")
bcd

#nmds code
set.seed(765)
nmds2 = metaMDS(m_com, distance = "bray")
nmds2
scores(nmds2)
en = envfit(nmds2, env, permutations = 999, na.rm = TRUE)
en
plot(nmds2)
plot(en)

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores = as.data.frame(scores(nmds2)$sites)
data.scores$RIVER = df$RIVER #adds in column with River 
data.scores$HABITAT = df$HABITAT #adds in column with Treatment
head(data.scores)
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)

#plot
gg2 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(shape = HABITAT, colour = RIVER), size = 5, alpha = 0.85) +  
  scale_shape_manual(values = c(16, 17)) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size = 0.75, alpha = 0.5, colour = "grey30",
               arrow = arrow(length=unit(0.25,"cm"), type = "closed")) +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            nudge_x = -0.03, nudge_y = -0.02,size=3, label = row.names(en_coord_cont)) + 
  theme(axis.title = element_text(size = 10, colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30")) +
  theme(legend.key = element_rect(colour = "transparent", fill = "white"), legend.position =  "none",
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
        legend.direction = "vertical", legend.box = "horizontal",
        legend.title = element_text(size=10),
        legend.text=element_text(size=8)) +
  labs(colour = "River", shape = "Habitat")
gg3 <- gg2 +  scale_colour_manual(values = brewer.pal(4, "Paired")[c(2,4)])
gg3

#rbcL - taxonomy-free
df<-read.csv("Diatom_taxfree_rbcl_20221028.csv")
df
colnames(df)
en<-read.csv("Env_20221027.csv")
en
colnames(en)

df$HABITAT <- as.character(df$HABITAT)
df$HABITAT <- factor(df$HABITAT, levels=unique(df$HABITAT))
df$RIVER <- as.character(df$RIVER)
df$RIVER <- factor(df$RIVER, levels=unique(df$RIVER))
com = df[,5:221]
env = en[,4:12]
#convert com to a matrix
m_com = as.matrix(com)
bcd = vegdist(m_com, "bray")
bcd

#nmds code
set.seed(765)
nmds2 = metaMDS(m_com, distance = "bray")
nmds2
scores(nmds2)
en = envfit(nmds2, env, permutations = 999, na.rm = TRUE)
en
plot(nmds2)
plot(en)

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores = as.data.frame(scores(nmds2)$sites)
data.scores$RIVER = df$RIVER #adds in column with River 
data.scores$HABITAT = df$HABITAT #adds in column with Treatment
head(data.scores)
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)


gg4 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(shape = HABITAT, colour = RIVER), size = 5, alpha = 0.85) +  
  scale_shape_manual(values = c(16, 17)) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size = 0.75, alpha = 0.5, colour = "grey30",
               arrow = arrow(length=unit(0.25,"cm"), type = "closed")) +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            nudge_x = -0.02, nudge_y = -0.05,size=3, label = row.names(en_coord_cont)) + 
  theme(axis.title = element_text(size = 10, colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30")) 
gg5 <- gg4 +  scale_colour_manual(values = brewer.pal(4, "Paired")[c(2,4)]) +
  theme(legend.key = element_rect(colour = "transparent", fill = "white"), legend.position =  "none",
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
        legend.direction = "vertical", legend.box = "horizontal",
        legend.title = element_text(size=10),
        legend.text=element_text(size=8)) +
  labs(colour = "River", shape = "Habitat")
gg5

#18S - taxonomy based
df<-read.csv("Diatom_taxbased_18S.csv")
df
colnames(df)
en<-read.csv("Env_20221027.csv")
en
colnames(en)

df$HABITAT <- as.character(df$HABITAT)
df$HABITAT <- factor(df$HABITAT, levels=unique(df$HABITAT))
df$RIVER <- as.character(df$RIVER)
df$RIVER <- factor(df$RIVER, levels=unique(df$RIVER))
com = df[,5:20]
env = en[,4:12]
#convert com to a matrix
m_com = as.matrix(com)
bcd = vegdist(m_com, "bray")
bcd

#nmds code
set.seed(765)
nmds2 = metaMDS(m_com, distance = "bray")
nmds2
scores(nmds2)
en = envfit(nmds2, env, permutations = 999, na.rm = TRUE)
en
plot(nmds2)
plot(en)

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores = as.data.frame(scores(nmds2)$sites)
data.scores$RIVER = df$RIVER #adds in column with River 
data.scores$HABITAT = df$HABITAT #adds in column with Treatment
head(data.scores)
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)

#plot
gg7 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(shape = HABITAT, colour = RIVER), size = 5, alpha = 0.85) +  
  scale_shape_manual(values = c(16, 17)) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size = 0.75, alpha = 0.5, colour = "grey30",
               arrow = arrow(length=unit(0.25,"cm"), type = "closed")) +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            nudge_x = 0.05, nudge_y = -0.02,size=3, label = row.names(en_coord_cont)) + 
  theme(axis.title = element_text(size = 10, colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30")) +
  theme(legend.key = element_rect(colour = "transparent", fill = "white"), legend.position = "none",
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
        legend.direction = "vertical", legend.box = "horizontal",
        legend.title = element_text(size=10),
        legend.text=element_text(size=8)) +
  labs(colour = "River", shape = "Habitat")
gg8 <- gg7 +  scale_colour_manual(values = brewer.pal(4, "Paired")[c(2,4)])
gg8

#18S - taxonomy-free
df<-read.csv("Diatom_taxfree_18S_20221028.csv")
df
colnames(df)
en<-read.csv("Env_20221027.csv")
en
colnames(en)

df$HABITAT <- as.character(df$HABITAT)
df$HABITAT <- factor(df$HABITAT, levels=unique(df$HABITAT))
df$RIVER <- as.character(df$RIVER)
df$RIVER <- factor(df$RIVER, levels=unique(df$RIVER))
com = df[,5:156]
env = en[,4:12]
#convert com to a matrix
m_com = as.matrix(com)
bcd = vegdist(m_com, "bray")
bcd

#nmds code
set.seed(765)
nmds2 = metaMDS(m_com, distance = "bray")
nmds2
scores(nmds2)
en = envfit(nmds2, env, permutations = 999, na.rm = TRUE)
en
plot(nmds2)
plot(en)

#extract NMDS scores (x and y coordinates) for sites from newer versions of vegan package
data.scores = as.data.frame(scores(nmds2)$sites)
data.scores$RIVER = df$RIVER #adds in column with River 
data.scores$HABITAT = df$HABITAT #adds in column with Treatment
head(data.scores)
en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)


gg9 = ggplot(data = data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(data = data.scores, aes(shape = HABITAT, colour = RIVER), size = 5, alpha = 0.85) +  
  scale_shape_manual(values = c(16, 17)) +
  geom_segment(aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2), 
               data = en_coord_cont, size = 0.75, alpha = 0.5, colour = "grey30",
               arrow = arrow(length=unit(0.25,"cm"), type = "closed")) +
  geom_text(data = en_coord_cont, aes(x = NMDS1, y = NMDS2), colour = "grey30", 
            nudge_x = -0.04, nudge_y = -0.065,size=3, label = row.names(en_coord_cont)) + 
  theme(axis.title = element_text(size = 10, colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30")) 
gg10 <- gg9 +  scale_colour_manual(values = brewer.pal(4, "Paired")[c(2,4)]) +
  theme(legend.key = element_rect(colour = "transparent", fill = "white"), legend.position =  "none",
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
        legend.direction = "vertical", legend.box = "horizontal",
        legend.title = element_text(size=10),
        legend.text=element_text(size=8)) +
  labs(colour = "River", shape = "Habitat")
gg10


