##EXTRACTING COLUMN NAMES (ALL EC #S)
library(stringr)
library(dplyr)
library(janitor)
library(tidyverse)
library(devtools); install_github('vqv/ggbiplot')
library(ggbiplot)
##CREATES HEADERS
bacteria_ec <- read.delim("bacteria_big_matrix_tab.txt", header = TRUE, quote = "",  sep= '\t')

header<-names(bacteria_ec)
complete_header <- gsub("X", "", header)
ec_names_only <- complete_header[c(-1)]

# # Read in the Processed Distance Matrix
# CLASS --> creates the dendrogram using the euclidean distance matrix
d_class <-read.delim('DistanceListforSynbio_Class.txt',header=FALSE, sep = "\t", row.names = 1)
colnames(d_class)<-c(complete_header)
names(d_class) <- row.names(d_class)
# Add zeroes along the diag back in
d_class[d_class<0.001] <- 0
d_class_asdistance <- as.dist(d_class)
fit_class<- hclust(d_class_asdistance, method="ward")
pdf("20220601 Class Tree.pdf", height = 8, width = 25)
plot(fit_class)
dev.off()


## COMPLETES PCA AND PLOTS PCA SPACE USING THE EC PRESENCE MATRIX 
bacteria_class <-read.delim("ec_space_class.txt", header = FALSE, quote = "",  sep= ' ')
colnames(bacteria_class)<-c(ec_names_only)
list_of_class<-row.names(d_class)
row.names(bacteria_class)<-list_of_class
ec_numbers <- ncol(bacteria_class)
bacteria_class_sum <- colSums(bacteria_class%>%select_if(is.numeric))
bacteria_class_sum[bacteria_class_sum>400]
sum(bacteria_class_sum/25899 > 0.75)
sum(bacteria_class_sum/25899 > 0.50)
sum(bacteria_class_sum/25899 > 0.25)
fit_prin_class <- princomp(t(bacteria_class[bacteria_class_sum>100]), cor=TRUE)  #data type is list
summary(fit_prin_class) # print variance accounted for
loadings(fit_prin_class) # pc loadings
#plot(fit_prin_class$scores[,2], fit_prin_class$scores[,3],type="lines") # scree plot
plot(fit_prin_class,type="lines") # scree plot
write.table(fit_prin_class$scores,"class_scores.txt") # the principal components
##NEED TO CHANGE THE COMPONENTS BEING PLOTTED
ggbiplot(fit_prin_class,choices=2:3, labels = row.names(fit_prin_class$scores))  #how to change the components tested
#biplot(lapply(fit_prin_class, "[", 2), lapply(fit_prin_class, "[", 3))


##CREATES HISTOGRAMS BASED ON THE TOP 20 EC CATEGORIES AND BOTTOM 20 
sorted_ec <- as.data.frame(cbind(sort(colSums(bacteria_class), T))) 
sorted_ec_noindex <- as.data.frame(cbind(row.names(sorted_ec),sorted_ec))
colnames(sorted_ec_noindex)<-c("EC Numbers", "EC Frequency")
row_sub = apply(sorted_ec, 1 , function(row) all(row!=0))
non_zerosorted <- as.data.frame(sorted_ec[row_sub,])
colnames(non_zerosorted)<- c("EC Frequency")
sorted_nonzero <-unique(merge(sorted_ec_noindex, non_zerosorted, by="EC Frequency"))
df<-subset(sorted_nonzero, select=c("EC Numbers", "EC Frequency"))
top_20 <- df %>% top_n(20)
barplot(top_20$`EC Frequency`, main = "Frequency of Top 20 ECs", xlab="EC Numbers" ,names.arg = top_20$`EC Numbers`)
bottom_20 <- df %>% top_n(-100)
unique_bottom <- bottom_20[!duplicated(bottom_20[ ,c("EC Frequency")]),]
barplot(unique_bottom$`EC Frequency`, main="Frequency of Bottom 23 ECs", xlab='EC Numbers', names.arg = unique_bottom$`EC Numbers`)
