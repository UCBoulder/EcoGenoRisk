library(dplyr)
library(ggplot2)

setwd("C:/Users/ulano/OneDrive/Desktop/CMBM-GEM/R studio/CompetitorFind/DIFFERENT")
# Summary lists that contain individual scores, functional distance, file name, and type of pairs
scores<-list('Scores')
distances<-list('Distances')
type<-list('Distance_Type')
file<-list('From which file')

# Assigns a CompetitorFind score based on overlapping substrates for competition and mutualism 
score<-function(competition, mutualism1, mutualism2){
  score<-(competition-(mutualism1+mutualism2)/2)/100
  return(score)
}

#======================== Enzymatically Different =============================#
# Extracting the score and distance from enzymatically different genome pairs
# The steps taken for each genome pair type is the same
# Opens the list of summary files created in CompetitorFind_test_cases.py
# Each file contains a list of number of total overlap substrates, number of mutualistic substrates, 
# functional distance between the two genomes, and substrates found in modified pathway
list_of_pairs_analyzed<-read.csv('diffferent.txt',sep='\t')
# Iterates through each of the files
for (i in 1:nrow(list_of_pairs_analyzed)){
  file_name<-list_of_pairs_analyzed[i,1]
  # Reads the summary file for the genome pair
  score_matrix<-read.csv(file_name,sep='\t', header=TRUE)
  # Extracts the functional distance between two genomes
  indv_dist<-as.double(score_matrix[2,8])
  # Extracts number of total overlap substrates
  competition<-as.double(score_matrix[2,2])
  # Extracts number of mutualistic substrates
  mutualism1<-as.double(score_matrix[2,5])
  mutualism2<-as.double(score_matrix[2,6])
  # Sends the substrate data for scoring
  indv_score<-score(competition, mutualism1, mutualism2)
  # Iteratively elongates all of the summary lists used for plotting
  # Appends the individual score to the overall score list
  scores<-append(scores, indv_score)
  # Appends individual distance to the overall distance list
  distances<- append(distances,indv_dist)
  # Appends the file name to the list of file names
  file<-append(file,file_name)
  # Appends the type of pairs analyzed to the overall types list
  type<-append(type,'different')
}

#======================== Enzymatically Similar =============================#
# Extracting the score and distance from enzymatically different genome pairs

#Changes directory where the similar files are saved
setwd("C:/Users/ulano/OneDrive/Desktop/CMBM-GEM/R studio/CompetitorFind/SIMILAR")
# Opens the list of summary files created in CompetitorFind_test_cases.py
# Each file contains a list of number of total overlap substrates, number of mutualistic substrates, 
# functional distance between the two genomes, and substrates found in modified pathway
list_of_pairs_analyzed<-read.csv('similar.txt',sep='\t')
# Iterates through each of the files
for (i in 1:nrow(list_of_pairs_analyzed)){
  file_name<-list_of_pairs_analyzed[i,1]
  # Reads the summary file for the genome pair
  score_matrix<-read.csv(file_name,sep='\t', header=TRUE)
  # Extracts the functional distance between two genomes
  indv_dist<-as.double(score_matrix[2,8])
  # Extracts number of total overlap substrates
  competition<-as.double(score_matrix[2,2])
  # Extracts number of mutualistic substrates
  mutualism1<-as.double(score_matrix[2,5])
  mutualism2<-as.double(score_matrix[2,6])
  # Sends the substrate data for scoring
  indv_score<-score(competition, mutualism1, mutualism2)
  # Iteratively elongates all of the summary lists used for plotting
  # Appends the individual score to the overall score list
  scores<-append(scores, indv_score)
  # Appends individual distance to the overall distance list
  distances<- append(distances,indv_dist)
  # Appends the file name to the list of file names
  file<-append(file,file_name)
  # Appends the type of pairs analyzed to the overall types list
  type<-append(type,'similar')
}


#======================== Enzymatically Random =============================#
# Extracting the score and distance from enzymatically different genome pairs

#Changes directory where the similar files are saved
setwd("C:/Users/ulano/OneDrive/Desktop/CMBM-GEM/R studio/CompetitorFind/RANDOM")
# Opens the list of summary files created in CompetitorFind_test_cases.py
# Each file contains a list of number of total overlap substrates, number of mutualistic substrates, 
# functional distance between the two genomes, and substrates found in modified pathway
list_of_pairs_analyzed<-read.csv('random.txt',sep='\t')
# Iterates through each of the files
for (i in 1:nrow(list_of_pairs_analyzed)){
  file_name<-list_of_pairs_analyzed[i,1]
  # Reads the summary file for the genome pair
  score_matrix<-read.csv(file_name,sep='\t', header=TRUE)
  # Extracts the functional distance between two genomes
  indv_dist<-as.double(score_matrix[2,8])
  # Extracts number of total overlap substrates
  competition<-as.double(score_matrix[2,2])
  # Extracts number of mutualistic substrates
  mutualism1<-as.double(score_matrix[2,5])
  mutualism2<-as.double(score_matrix[2,6])
  # Sends the substrate data for scoring
  indv_score<-score(competition, mutualism1, mutualism2)
  # Iteratively elongates all of the summary lists used for plotting
  # Appends the individual score to the overall score list
  scores<-append(scores, indv_score)
  # Appends individual distance to the overall distance list
  distances<- append(distances,indv_dist)
  # Appends the file name to the list of file names
  file<-append(file,file_name)
  # Appends the type of pairs analyzed to the overall types list
  type<-append(type,'random')
}

#========================= Dataframe Formatting   =============================#

# Binds the summary lists together, forming a list with four columns
new<-cbind(unlist(type), unlist(scores), unlist(distances), unlist(file))
# Removes the initial row used for initializing the lists
scoring_results<-new[-1,]
# Forms the four column list into a dataframe
scoring_results<-as.data.frame(new)
# Sets appropriate column names
colnames(scoring_results)<-c('Type','Scores','Distance','Filename')
# Sets the type column as a factor since the datapoints will be colored based on type
scoring_results$Type<-factor(scoring_results$Type)
# Sets scores and distance as numeric types 
scoring_results$Scores<-as.numeric(scoring_results$Scores)
scoring_results$Distance<-as.numeric(scoring_results$Distance)
# Removes any rows that got misread during data extraction
scoring_results<-scoring_results[! grepl('Name_of_Genome', scoring_results$Distance),]
scoring_results<-scoring_results[! grepl('Name_of_Genome', scoring_results$Scores),]


# Convert Scores and Distance to numeric and remove rows with missing values
scoring_results <- scoring_results %>%
  mutate(
    Scores = as.numeric(Scores),
    Distance = as.numeric(Distance)
  ) %>%
  na.omit()

#=========================== Plotting Results   ===============================#

# Create a ggplot object with a scatter plot layer
p <- ggplot(scoring_results, aes(x = Distance, y = Scores, color = Type)) +
  geom_point()

# Set the x and y axis limits and tick marks
p <- p + scale_x_continuous(limits = c(0, 24), breaks = seq(0, 24, 2)) +
  scale_y_continuous(limits = c(0, 5), breaks = seq(0, 5, 0.5))

# Set the labels for the axes and legend
p <- p + labs(x = "Distance", y = "Scores", color = "Type")


# Set the colors for the legend
p <- p + scale_color_manual(values = c("different" = "red", "similar" = "blue", "random" = "chartreuse4"))

# Display the plot
p
