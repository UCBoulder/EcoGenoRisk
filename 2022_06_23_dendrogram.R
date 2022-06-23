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

lineage <-read.delim('cleaned_bac120_taxonomy_r207.tsv.txt',header=TRUE, sep = "\t")
bacteria_matrix <-  read.table("bacteria_combined1.csv", header=TRUE, quote = "", sep = ' ')
lineage_isolated <-data.frame(lineage$GCF, lineage$Class, lineage$Phylum)
colnames(lineage_isolated)<- c("GCF", "Class", "Phylum")
bacteria_isolated<-data.frame(bacteria_matrix$Name_of_Genome)
colnames(bacteria_isolated)<-c("GCF")
bacteria_split<- str_split_fixed(bacteria_isolated$GCF, '_', 3)
bacteria_isolated$GCF<-paste("GCF_",bacteria_split[,2], sep="")
merged_lineage<- merge(bacteria_isolated, lineage_isolated, by="GCF", all.x=TRUE)
reference_key<-data.frame(merged_lineage$Class, merged_lineage$Phylum)
colnames(reference_key)<-c("Class", "Phylum")

d_class <-read.delim('no_filter_DistanceListforSynbio_Class.txt',header=FALSE, sep = "\t", row.names = 1)
names(d_class) <- row.names(d_class)
#reference_key <- read.delim("class_phylum_key.txt", header=TRUE, sep='\t')

##DENDROGRAMS AND PHYLA COLORING: https://stackoverflow.com/questions/43145448/how-to-color-dendrogram-labels-using-r-based-on-label-name-not-grouping
d_class_factor<-factor(reference_key$Phylum)
d_class[(d_class)<0.001] <- 0
d_class_asdistance <- as.dist(d_class)
fit_class<- hclust(d_class_asdistance, method="ward.D2")
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