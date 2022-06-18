setwd("C:/Users/ulano/OneDrive/Desktop/CMBM-GEM/R studio")
install.packages(c("ggplot2", "ggdendro", "janitor", "devtools","tidyverse", "stringr", "ggraph","dendextend", "randomcoloR"))
install.packages("dbplyr", type = "binary")
library(stringr)
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
dend1 <- as.dendrogram(fit_class)
dend2 <- as.dendrogram(fit_class)
dend3 <- as.dendrogram(fit_class)

##METHOD 1:  https://stackoverflow.com/questions/27485549/how-to-colour-the-labels-of-a-dendrogram-by-an-additional-factor-variable-in-r
colors_to_use <- as.numeric(d_class_factor)
#colors_to_use <- colors_to_use[order.dendrogram(dend1)]
labels_colors(dend1)<-colors_to_use
labels_colors(dend1)
color_map <-as.data.frame(labels_colors(dend1))
write.table(color_map, "color_map.txt", sep='\t')
color_map<-as.data.frame(labels_colors(dend1))
plot(dend1, main="Method 1")
par(mar = c(9, 0.5, 0.5, 0.5)) # Set the margin on all sides to 6
legend(x="bottom", inset=c(0,-1),legend=unique(d_class_factor), col =unique(colors_to_use))


##METHOD 2:  https://stackoverflow.com/questions/41704404/r-plot-color-legend-by-factor
n <- 42
palette <- distinctColorPalette(n)
for (i in 1:42){
  rbg_conversion[i]<-paste("c(",paste(as.vector(col2rgb(palette[i])), collapse = ","), ")")
}
labels_colors(dend2)<-palette
par(mar = c(9, 0.5, 0.5, 0.5)) # Set the margin on all sides to 6
dend1 <- as.dendrogram(fit_class)
plot(dend2, pch=16, col=d_class_factor[["Phylum"]], main="Method 2")
legend("bottom", inset = c(0, -0.5), legend=unique(d_class_factor[["Phylum"]]), col=rbg_conversion,  xpd = TRUE, horiz = TRUE)
#legend("topleft", legend=unique(d_class_factor), pch=16, col=unique(d_class_factor))

##METHOD 3: https://stackoverflow.com/questions/8045538/labelling-ggdendro-leaves-in-multiple-colors
d_class_factor3<-factor(reference_key)
p2<- ggplot(segment(dendro_data(dend3)))+geom_segment(aes(x=d_class_asdistance, y=0))
p3<-p2+geom_text(data=label(dendro_data(dend3)), aes(label = label, x=d_class_asdistance, y=0, col=d_class_factor3[["Phylum"]]))
plot(p3, main="Method 3")

