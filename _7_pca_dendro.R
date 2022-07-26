#setwd("C:/Users/ulano/OneDrive/Desktop/CMBM-GEM/R studio/Documents Required for PCA and Dendrogram Runs/Fungi")
# If there are issues with installing/reinstalling packages, then use the following link and code below as a guide: 
# https://stackoverflow.com/questions/63390194/package-xxx-was-installed-before-r-4-0-0-please-re-install-it
##================UNINSTALLING AND REINSTALLING PACKAGES======================##
# # check your package library path
# .libPaths()
# 
# # grab old packages names
# old_packages <- installed.packages(lib.loc = "/Library/Frameworks/R.framework/Versions/3.6/Resources/library")
# old_packages <- as.data.frame(old_packages)
# list.of.packages <- unlist(old_packages$Package)
# 
# # remove old packages
# remove.packages( installed.packages( priority = "NA" )[,1] )
# 
# # reinstall all packages
# new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
# if(length(new.packages)) install.packages(new.packages)
# lapply(list.of.packages,function(x){library(x,character.only=TRUE)})
##============================================================================##
update.packages(checkBuilt = TRUE, ask = FALSE)

setwd("C:/Users/ulano/OneDrive/Desktop/CMBM-GEM/R studio/Documents Required for PCA and Dendrogram Runs/Fungi")
install.packages(c("ggplot2", "ggdendro", "janitor", "devtools","tidyverse", "stringr", "ggraph","dendextend", "randomcoloR", "colorspace"))
install.packages("dbplyr", type = "binary")
install.packages("ade4") 
library(stringr)
library(colorspace)
library(janitor)
library(tidyverse)
library(ggbiplot)
library(ggdendro)
library(dendextend)
library(dbplyr)
library(randomcoloR)
library(ade4)

# name1<-'taxonomy.tsv'
# name2<-"complete_binary_matrix.txt"
# name3<-'ec_space_class.txt'
# name4<-'Distance_matrix_Combined.txt'
# name5 <-'wfilter_Distance_matrix_Combined.txt'

rendering_plots<-function(name1, name2, name3, name4, loc){
  ##=========================READING IN FILES===================================##
  # Output from python  script 0, contains the full lineage for both archaea
  # and bacteria, along with GCF annotation
  lineage <-read.delim(name3, header=TRUE, sep = "\t")
  
  #Output from python script 0, contains the binary matrices for archaea and bacteria
  binary_matrix <-  read.table(paste(loc, name1, sep = '/'), header=TRUE, row.names = 'Name_of_Genome', quote = "", sep = '\t')
  
  # Output from python script 5, contains the EC space after binning depending on which rank the user selected
  ec_rank <-read.delim(paste(loc, name2, sep = '/'), header = TRUE ,quote = "",  sep= '\t')
  # Distance matrix from script 5, used for dendrogram construction 
  d_rank <-read.delim(paste(loc,name4, sep = '/'), header=FALSE, quote = "", sep = "\t")
  print('Document upload is complete')
  ##=====================FINDS RANK USED FOR CLUSTERING=========================##
  # Uses EC averaged space for finding the rank
  rank_txt <- unlist(strsplit(name2,'_'))[3]
  rank <- str_to_title(unlist(strsplit(rank_txt, '.txt')))
  index<-grep(rank, names(lineage), ignore.case ="True")
  rank_above <-names(lineage)[index-1]
  ##=========================DATAFRAME FORMATTING===============================##
  # Checks that there are no row repeat in the overall combined matrix, can cause
  # issues in future processing. If repeats are found, then those rows are dropped
  # Returns the list of genome GCFs that are repeated
  n_occur <- data.frame(table(binary_matrix$Name_of_Genome))
  n_occur[n_occur$Freq > 1,]
  repeated_rows<-binary_matrix[binary_matrix$Name_of_Genom %in% n_occur$Var1[n_occur$Freq > 1],]
  binary_matrix<-binary_matrix[!(binary_matrix$Name_of_Genome %in% c(repeated_rows)), ]
  
  # Creates a reference dataframe for dendrogram coloring-- not complete
  # Creates a df containing only the GCF annotation, class and phylum 
  reference_key <-data.frame(lineage[,rank], lineage[,rank_above])
  # Names the columns for consistency
  colnames(reference_key)<-c(rank, rank_above)
  
  row.names(ec_rank)<-ec_rank[,rank]
  ec_rank[,rank]<-NULL
  # Removes X in front of the EC numbers
  colnames(ec_rank) <- gsub("X", "", names(ec_rank))
  
  # Sets the square distance matrix column and index names as the rank names
  row.names(d_rank)<-row.names(ec_rank)
  colnames(d_rank)<-row.names(ec_rank)
  
  ##=====================COMPLETES AND PLOTS PCA================================##
  # Uses the phylogenetically clustered EC matrix, removes data that provides no 
  # variability
  # Formats the documents to have correct headers and indexes 
  # Finds the names of the ranks, lineage bin names 
  # Number of ranks present, how many bins were created during clustering in script 5
  n<-nrow(ec_rank)
  nonzero_ec_num <- ec_rank%>% select_if(colSums(ec_rank)!=0) 
  summations_for_ec <- colSums(nonzero_ec_num)
  
  #summations_for_ec<-colSums(ec_rank)
  # Finds  three bins where sums are greater than 400
  sum(summations_for_ec/n > 0.75)
  sum(summations_for_ec/n > 0.50)
  sum(summations_for_ec/n > 0.25)
  
  # Creates a principle component table
  fit_prin_rank <- princomp(ec_rank[colSums(ec_rank)/n > 0.25], cor= TRUE)  #data type is list
  summary(fit_prin_rank) # print variance accounted for
  loadings(fit_prin_rank) # pc loadings
  plot(fit_prin_rank,type="lines") # scree plot
  write.table(fit_prin_rank$scores,"rank_scores.txt") # the principal components
  
  # Plots PCA table by graphing component 2 and 3
  # If encountering issues with ggbiplot invalid rot value, then follow this guide: 
  # https://stackoverflow.com/questions/27016619/prcomp-and-ggbiplot-invalid-rot-value
  ggbiplot(fit_prin_rank,choices = 2:3, labels = row.names(fit_prin_rank$scores))  
  
  ##===========SUMMARY TABLES OF TOP 20 EC CATEGORIES AND BOTTOM 20=============##
  ## Uses unbinned and unweighted EC space
  #Finds the total number of samples used
  sample_number <- nrow(binary_matrix)
  
  #Finds column sums for each EC number, divides by the total number of organisms to find percent frequency
  # creates a dataframe of EC number sorted sums
  binary_matrix<-select_if(binary_matrix, is_numeric)
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
  
  ##======================DENDROGRAMS AND PHYLA COLORING========================##
  # Uses the distance matrix post phylogenetic clustering, colors data based on the
  # Rank above the one chosen 
  # Use link below as a guide for trouble shooting
  # https://stackoverflow.com/questions/43145448/how-to-color-dendrogram-labels-using-r-based-on-label-name-not-grouping
  #Creates a dendrogram factor
  d_rank_factor<-factor(reference_key[,rank_above])
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
  plot(dend4)
  #dev.off()
  #Creates a distinct color palette for each Phyla
  #palette <- distinctColorPalette(n)
  #cols<-rainbow_hcl(number_group)
  #col_group<- cols[(d_rank_factor)]
  # Creates a matrix which assigns a color for each Phylum
  #col_group <-palette[d_rank_factor[['Phylum']]]
  #col_group<- col_group[order.dendrogram(dend4)]
  #Sets margin size
  par(mar = c(20, 1, 1, 1))
  #Customizes dendrogram by coloring groups, setting label size, and makes the labels uniform as the individual branch hieght gets adjusted
  #dend4 <- dend4%>% set("labels_colors", col_group[['Phylum']])%>% set("labels_cex", 1) %>% raise.dendrogram (-1) %>% plot(main ="Method") # Set the margin on all sides to 6
  #legend("bottom", legend = levels(unique(d_rank_factor[['Phylum']])),  fill=palette, cex=0.6, x.intersp =0.1, y.intersp=0.75, title= "Associated Phylum", inset= c(0,0), ncol=3) 
  
  ##============================= MANTEL TEST ==================================##
  weighted_distances<- read.delim('wfilter_Distance_matrix_Combined.txt',header=FALSE, quote = "", sep = "\t")
  unweighted_distances <- read.delim('Distance_matrix_Combined.txt', header=FALSE, quote = "", sep = "\t")
  weighted_distances_dist<-dist(weighted_distances)
  unweighted_distances_dist<-dist(unweighted_distances)
  mantel.rtest(m1=weighted_distances_dist,m2=unweighted_distances_dist, nrepet = 9999)
  
  
}

