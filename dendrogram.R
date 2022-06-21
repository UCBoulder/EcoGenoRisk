setwd("C:/Users/ulano/OneDrive/Desktop/CMBM-GEM/R studio")
install.packages(c("ggplot2", "ggdendro", "janitor", "devtools","tidyverse", "stringr", "ggraph","dendextend", "randomcoloR", "colorspace"))
install.packages("dbplyr", type = "binary")
library(stringr)
library(colorspace)
library(janitor)
library(tidyverse)
library(devtools); install_github('vqv/ggbiplot')
library(ggbiplot)
library(ggdendro)
library(dendextend)
library(dbplyr)
library(randomcoloR)

d_class <-read.delim('no_filter_DistanceListforSynbio_Class.txt',header=FALSE, sep = "\t", row.names = 1)
names(d_class) <- row.names(d_class)

reference_key <- read.delim("class_phylum_key.txt", header=TRUE, sep='\t')
d_class_factor<-factor(reference_key$Phylum)
d_class[(d_class)<0.001] <- 0
d_class_asdistance <- as.dist(d_class)
fit_class<- hclust(d_class_asdistance, method="ward.D2")

##METHOD 4: https://stackoverflow.com/questions/43145448/how-to-color-dendrogram-labels-using-r-based-on-label-name-not-grouping
dend4 <- as.dendrogram(fit_class)
number_group<-length(unique(d_class_factor))
n <- number_group
palette <- distinctColorPalette(n)
#cols<-rainbow_hcl(number_group)
#col_group<- cols[(d_class_factor)]
col_group <-palette[d_class_factor]
col_group<- col_group[order.dendrogram(dend4)]
par(mar = c(20, 1, 1, 1))
dend4 <- dend4%>% set("labels_colors", col_group)%>% set("labels_cex", 1) %>% raise.dendrogram (-1) %>% plot(main ="Method") # Set the margin on all sides to 6
legend("bottom", legend = levels(unique(d_class_factor)),  fill=palette, cex=0.6, x.intersp =0.1, y.intersp=0.75, title= "Associated Phylum", inset= c(0,0), ncol=3) 
#inset= c(1,1)