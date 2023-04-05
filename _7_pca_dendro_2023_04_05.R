# If there are issues with installing/reinstalling packages, then use the following link and code below as a guide:
# https://stackoverflow.com/questions/63390194/package-xxx-was-installed-before-r-4-0-0-please-re-install-it
#================UNINSTALLING AND REINSTALLING PACKAGES======================##
# check your package library path
.libPaths()

# grab old packages names
old_packages <- installed.packages(lib.loc = "C:/Users/ulano/AppData/Local/R/win-library/4.2")
old_packages <- as.data.frame(old_packages)
list.of.packages <- unlist(old_packages$Package)

# remove old packages
remove.packages( installed.packages( priority = "NA" )[,1] )

# reinstall all packages
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
lapply(list.of.packages,function(x){library(x,character.only=TRUE)})
##============================================================================##
#update.packages(checkBuilt = TRUE, ask = FALSE)

setwd("C:/Users/ulano/OneDrive/Desktop/CMBM-GEM/R studio/2023_3_17_Results/Results_Docs")
install.packages(c("ggplot2", "ggdendro", "janitor", "usethis", "devtools", "dplyr", 'ggplot2'
                   ,"tidyverse", "circlize", "stringr", "ggraph","dendextend", "randomcoloR", "colorspace", dependencies = TRUE))
install.packages("dbplyr", type = "binary", dependencies = TRUE)
install.packages("ade4")
install.packages("tools")
library(usethis)
library(devtools)
#install_github("vqv/ggbiplot")
library(tools); packageVersion('tools')
library(stringr); packageVersion('stringr') #V 1.4.0
library(colorspace); packageVersion('colorspace') #V 2.0.3
library(janitor); packageVersion('janitor') #
library(tidyverse); packageVersion('tidyverse') #V
library(ggbiplot); packageVersion('ggbiplot')
library(ggdendro); packageVersion('ggdendro')
library(dendextend); packageVersion('dendextend')
library(dbplyr); packageVersion('dbplyr')
library(randomcoloR); packageVersion('randomcoloR')
library(ade4); packageVersion('ade4') #V 1.7.19
library(circlize); packageVersion('circlize')
library(dplyr); packageVersion('dplyr') #V 1.7.19
library(ggplot2); packageVersion('ggplot2')

name3<-'taxonomy_2023_1_5.tsv'
name1<-"combined_binary_summary_matrix_2023_3_15.csv"
name2<-'dendro_ec_space_clustered_unweightedclass (1).txt'
name4<-'class_grouped_distance_matrix.txt'

##=========================READING IN FILES===================================##
# Output from python  script 0, contains the full lineage for both archaea
# and bacteria, along with GCF annotation
lineage <-read.delim(name3, header=TRUE, sep = "\t")

#Output from python script 0, contains the binary matrices for archaea and bacteria
binary_matrix <-  read.table(name1, header=TRUE, quote = "", sep = '\t')

# Output from python script 5, contains the EC space after binning depending on which rank the user selected
ec_rank <-read.delim(name2, header = TRUE ,quote = "",  sep= '\t', row.names = 1)
# Distance matrix from script 5, used for dendrogram construction 
d_rank <-read.delim(name4, header=TRUE, quote = "", sep = "\t", row.names = 1)

print('Document upload is complete')
##===================== RANK USED FOR CLUSTERING=========================##
rank<-'Class'
rank_above<-'Phylum'
##=========================DATAFRAME FORMATTING===============================##
# Checks that there are no row repeat in the overall combined matrix, can cause
# issues in future processing. If repeats are found, then those rows are dropped
# Returns the list of genome GCFs that are repeated
n_occur <- data.frame(table(binary_matrix$Name_of_Genome))
n_occur[n_occur$Freq > 1,]
repeated_rows<-binary_matrix[binary_matrix$Name_of_Genom %in% n_occur$Var1[n_occur$Freq > 1],]
binary_matrix<-binary_matrix[!(binary_matrix$Name_of_Genome %in% c(repeated_rows)), ]

# Removes all X from the column names
colnames(ec_rank) <- gsub("X", "", names(ec_rank))


##=====================COMPLETES AND PLOTS PCA================================##
# Uses the phylogenetically clustered EC matrix, removes data that provides no 
# variability
# Formats the documents to have correct headers and indexes 
# Finds the names of the ranks, lineage bin names 
# Number of ranks present, how many bins were created during clustering in script 5
n<-nrow(ec_rank)
nonzero_ec_num <- ec_rank%>% select_if(colSums(ec_rank)/n>0.01) 
summations_for_ec <- colSums(nonzero_ec_num)

# Transposes the nonzero EC BSM
transposed<-t(nonzero_ec_num)
# Completes principle component analysis
pca_transposed<-prcomp(transposed)

# Finds the variance explained by each component
percent_variance <- round(100 * (summary(pca_transposed)$importance[2,]), 2)
pc1_var <- paste0("PC2 (", percent_variance[2], "%)")
pc2_var <- paste0("PC3 (", percent_variance[3], "%)")

# Creates a biplot for PC2 vs. PC3
biplot(pca_transposed, choices=2:3, xlab = pc1_var, ylab = pc2_var,  col=c("black","chartreuse4"))

# Finds the variance explained for each component for Scree plot
var_explained=(pca_transposed$sdev)^2/sum((pca_transposed$sdev)^2)
var_explained_df <- data.frame(PC= paste0("PC",1:5), Variance_Explained=var_explained[1:5])
head(var_explained_df)
# Creates Scree Plot
var_explained_df %>%
  ggplot(aes(x=PC,y=Variance_Explained, group=1))+
  geom_point(size=4)+
  geom_line()+
  labs(title="Scree plot of Class-Grouped EC BSM PCA")+ 
  xlab('Components')+
  ylab('Variance Explained')

# Creates additional NMDS for analysis

library(vegan)

nmds<-metaMDS(nonzero_ec_num, k=2, trymax=1000)
stressplot(nmds)
plot(nmds)
ordiplot(nmds,type="n")
orditorp(nmds,display='species',col="red",air=0.01)
orditorp(nmds,display='sites',cex=1,air=0.01)


##===========SUMMARY TABLES OF TOP 20 EC CATEGORIES AND BOTTOM 20=============##
## Uses unbinned and unweighted EC space
#Finds the total number of samples used
sample_number <- nrow(binary_matrix)

#Finds column sums for each EC number, divides by the total number of organisms to find percent frequency
# creates a dataframe of EC number sorted sums
sorted_ec <- as.data.frame(cbind(sort(colSums(binary_matrix[,-1])/sample_number, T)))

# Binds the EC numbers to the sums
sorted_ec_noindex <- as.data.frame(cbind(row.names(sorted_ec),sorted_ec))

# Adds column names
colnames(sorted_ec_noindex)<-c("EC Numbers", "Percent EC Frequency")

#Finds non zero sums for EC
row_sub = apply(sorted_ec, 1 , function(row) all(row!=0))

# Creates a data frame of the nonzero sums
non_zerosorted <- as.data.frame(sorted_ec[row_sub,])
colnames(non_zerosorted)<- c("Percent EC Frequency")
sorted_nonzero <-unique(merge(sorted_ec_noindex, non_zerosorted, by="Percent EC Frequency"))
df<-subset(sorted_nonzero, select=c("EC Numbers", "Percent EC Frequency"))
# Finds top 20 most frequent EC numbers which correlate to the essential gene set
top_20 <- df %>% top_n(20)
bottom_20 <- df %>% top_n(-700)

#Finds unique bottom ec numbers (700 were chosen since the lower the frequency, the more repetitive the number)
unique_bottom <- bottom_20[!duplicated(bottom_20[ ,c("Percent EC Frequency")]),]

# Saves as tables
write.table(unique_bottom, "bottom_23.txt")
write.table(top_20, "top_20.txt")

##======================DENDROGRAMS AND PHYLA COLORING========================##
# Uses the distance matrix post phylogenetic clustering, colors data based on the
# Rank above the one chosen 
# Use link below as a guide for trouble shooting
# https://stackoverflow.com/questions/43145448/how-to-color-dendrogram-labels-using-r-based-on-label-name-not-grouping
#Creates a dendrogram factor
d_rank[(d_rank)<0.001] <- 0
# Creates a distance matrix by finding the distance between rows of matrix
d_rank_asdistance <- as.dist(d_rank)
# Completes cluster analysis  by using the Ward method
# This is used for tree construction
fit_rank<- hclust(d_rank_asdistance, method="ward.D2")
#Creates dendrogram of the distance matrix
dend4 <- as.dendrogram(fit_rank)
# Creates circular dendrogram 
circlize_dendrogram(dend4, labels_track_height = 0.4, dend_track_height = 0.3)

##============================= Robinson-Foulds Metric ==================================##
library(ape)

setwd("C:/Users/ulano/OneDrive/Desktop/CMBM-GEM/R studio/2023_3_17_Results/Results_Docs/Robinson-Foulds")#weighted_distances<- read.delim('/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/Chimera1/Only one chimera/diff_weighted_clustered_distance_matrix.txt',header=FALSE, quote = "", sep = "\t")
# Reads in weighted, taxonomy-grouped EC space
weighted_distances<- read.delim('weighted_distance_matrix_for_dendro.txt',header=TRUE, row.names=1, quote = "", sep = "\t")
# Reads in unweighted, taxonomy-grouped EC space
unweighted_distances <- read.delim('distance_matrix_for_dendro.txt',header=TRUE, row.names=1, quote = "", sep = "\t")
# Reformats both as distance matrices
mat_weighted <- as.dist(weighted_distances)
mat_unweighted <- as.dist(unweighted_distances)
# Preforms hierarchical cluster analysis for each distance matrix
hc_weighted<-hclust(mat_weighted)
hc_unweighted<-hclust(mat_unweighted)
# Converts into a tree object "phylo" 
weighted_tree<-as.phylo(hc_weighted)
unweighted_tree<-as.phylo(hc_unweighted)
# Returns the "number of differences" 
rf<-dist.topo(weighted_tree,unweighted_tree)