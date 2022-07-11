setwd("/home/anna/rstudio-2022.02 (1).3-492-amd64-debian")
install.packages(c("ggplot2", "ggdendro", "janitor", "usethis","tidyverse", "stringr", "devtools", "ggraph","dendextend", "randomcoloR", "colorspace"))
install.packages("dbplyr", type = "binary")
library(stringr)
library(colorspace)
library(janitor)
library(tidyverse)
library(ggbiplot)
library(ggdendro)
library(dendextend)
library(dbplyr)
library(randomcoloR)

data_analysis_outputs<- function(name1, name2, name3){
  #Output from python  script 0, contains the full lineage for both archaea and bacteria, along with GCF annotation
  lineage <-read.delim('taxonomy.tsv',header=TRUE, sep = "\t")
  #Output from python script 0, contains the binary matrices for archaea and bacteria
  binary_matrix <-  read.table(name1, header=TRUE, quote = "", sep = '\t')
  # Output from python script 5, contains the EC space after binning depending on which rank the user selected
  ec_rank <-read.delim(name2, header = FALSE, quote = "",  sep= '\t')
  # Distance matrix from script 5, used for dendrogram construction 
  d_rank <-read.delim(name3, header=FALSE, sep = "\t", row.names = 1)
  
  ## Creates a reference dataframe for dendrogram coloring-- not complete
  # Creates a df containing only the GCF annotation, class and phylum 
  reference_key <-data.frame(lineage$Class, lineage$Phylum)
  # Names the columns for consistency
  colnames(reference_key)<-c("Class", "Phylum")
  
  ## FORMATTING BY CREATING HEADERS/ROW NAMES
  header<-names(binary_matrix)
  # Removes X from the EC numbers 
  complete_header <- gsub("X", "", header)
  ec_names_only <- complete_header[c(-1)]
  
  ## COMPLETES PCA AND PLOTS PCA SPACE USING THE EC PRESENCE MATRIX 
  # Formatting the documents to have correct headers and indexes 
  colnames(ec_rank)<-c(ec_names_only)
  # Finds the names of the ranks, lineage bin names 
  list_of_rank<-row.names(d_rank)
  # Sets index names to the distance matrix
  row.names(ec_rank)<-list_of_rank
  # Number of ranks present, how many bins were created during clustering in script 5
  ec_numbers <- ncol(ec_rank)
  #Sums the columns of EC space
  ec_rank_sum <- colSums(ec_rank)
  # Finds  three bins where sums are greater than 400
  ec_rank_sum[ec_rank_sum>400]
  sum(ec_rank_sum/99 > 0.75)
  sum(ec_rank_sum/99 > 0.50)
  sum(ec_rank_sum/99 > 0.25)
  # Creates a principle component table
  fit_prin_rank <- princomp(t(ec_rank[ec_rank_sum]), cor=TRUE)  #data type is list
  summary(fit_prin_rank) # print variance accounted for
  loadings(fit_prin_rank) # pc loadings
  plot(fit_prin_rank,type="lines") # scree plot
  write.table(fit_prin_rank$scores,"rank_scores.txt") # the principal components
  # Plots PCA table by graphing component 2 and 3
  ggbiplot(fit_prin_rank,choices=2:3, labels = row.names(fit_prin_rank$scores))  
  
  ##CREATES HISTOGRAMS BASED ON THE TOP 20 EC CATEGORIES AND BOTTOM 20 
  ## Uses unbinned and unweighted EC space
  #Finds the total number of samples used
  sample_number <- nrow(binary_matrix)
  #Finds column sums for each EC number, divides by the total number of organisms to find percent frequency
  # creates a dataframe of EC number sorted sums
  sorted_ec <- as.data.frame(cbind(sort(colSums(binary_matrix)/sample_number, T)))
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
  
  
  ##DENDROGRAMS AND PHYLA COLORING: https://stackoverflow.com/questions/43145448/how-to-color-dendrogram-labels-using-r-based-on-label-name-not-grouping
  #Creates a dendrogram factor
  d_rank_factor<-factor(reference_key)
  
  d_rank[(d_rank)<0.001] <- 0
  # Creates a distance matrix by finding the distance between rows of matrix
  d_rank_asdistance <- as.dist(d_rank)
  # Completes cluster analysis  by using the Ward method
  # This is used for tree construction
  fit_rank<- hclust(d_rank_asdistance, method="ward.D2")
  #Creates dendrogram of the distance matrix
  dend4 <- as.dendrogram(fit_rank)
  #Finds number of unique Phylums present
  number_group<-length(unique(d_rank_factor))
  n <- number_group
  # Saves dendrogram 1 as black and white diagram 
  pdf(dend4, "monochroma")
  dev.off()
  #Creates a distinct color palette for each Phyla
  palette <- distinctColorPalette(n)
  #cols<-rainbow_hcl(number_group)
  #col_group<- cols[(d_rank_factor)]
  # Creates a matrix which assigns a color for each Phylum
  col_group <-palette[d_rank_factor[['Phylum']]]
  col_group<- col_group[order.dendrogram(dend4)]
  #Sets margin size
  par(mar = c(20, 1, 1, 1))
  #Customizes dendrogram by coloring groups, setting label size, and makes the labels uniform as the individual branch hieght gets adjusted
  dend4 <- dend4%>% set("labels_colors", col_group[['Phylum']])%>% set("labels_cex", 1) %>% raise.dendrogram (-1) %>% plot(main ="Method") # Set the margin on all sides to 6
  legend("bottom", legend = levels(unique(d_rank_factor[['Phylum']])),  fill=palette, cex=0.6, x.intersp =0.1, y.intersp=0.75, title= "Associated Phylum", inset= c(0,0), ncol=3) 
}
