install.packages('ggplot2')
install.packages('dplyr')
install.packages('tidyverse')
install.packages("rstatix")
install.packages("ggpubr")

library(ggplot2)
library(dplyr)
library(tidyverse)
library(rstatix)
library(ggpubr)

# Changes the directory where method1_cummulative_results.csv is saved
setwd("C:/Users/ulano/OneDrive/Desktop/CMBM-GEM/R studio/Friedman/Method 1")
# Reads the summary document for Method 1. Document contains the genome group 
# pair tested, ecological relationship, substrate overlap, and the growth medium tested
method_1_results <-
  read.csv(
    'method1_cummulative_results.csv',
    header = TRUE,
    quote = "",
    sep = ','
  )
# Isolates the substrate overlap score and the ecological relationship 
method1 <-
  data.frame(method_1_results$Score, method_1_results$Ecology)
colnames(method1) <- c('Score', 'Ecological_Relationship')

#==============================================================================#
# Makes a dot plot with box plot overlap of the dataset, labels the dot plots,
# reorders the x-axis with the levels command 
p <-
  ggplot(method1, aes(
    x = factor(
      Ecological_Relationship,
      level = c(
        'Parasitism',
        'Amensalism',
        'Competition',
        'Neutralism',
        'Commensalism',
        'Mutualism',
        'Inconclusive'
      )
    ),
    y = Score
  )) +
  geom_dotplot(binaxis = 'y',
               stackdir = 'center' ,
               binwidth = 0.2) + labs(title = 'Method 1: Ecological Relationships By Score', x =
                                        'Ecological Relationship', y = 'Numnber of Substrates Overlapped')

# Calculates the mean substrate overlap for each ecological relationship
means <- aggregate(Score ~ Ecological_Relationship, method1, mean)

# Calculates the standard deviation for each ecological relationship
standard_dev <-
  aggregate(Score ~ Ecological_Relationship, method1, sd)

# Adds Wilcoxon relationships with Competition as the reference
p1 <- p +  geom_boxplot() +
  geom_dotplot(binaxis = 'y',
               stackdir = 'center',
               binwidth = 0.2) +
  stat_compare_means(aes(label = paste("p = ", ..p.format..)), method = "wilcox.test", ref.group = "Competition")

# Calculates the number of pairs that have a mutualistic relationship 
number_mutualism <-
  sum(method1$Ecological_Relationship == 'Mutualism')
# Calculates the number of pairs that have a competition relationship 
number_competition <-
  sum(method1$Ecological_Relationship == 'Competition')
# Calculates the number of pairs that have a neutralistic relationship 
number_neutralism <-
  sum(method1$Ecological_Relationship == 'Neutralism')
# Calculates the number of pairs that have a amensalism relationship 
number_amensalism <-
  sum(method1$Ecological_Relationship == 'Amensalism')
# Calculates the number of pairs that have a inconclusive relationship 
number_inconclusive <-
  sum(method1$Ecological_Relationship == 'Inconclusive')
# Calculates the number of pairs that have a parasitism relationship 
number_parasitism <-
  sum(method1$Ecological_Relationship == 'Parasitism')
# Calculates the number of pairs that have a commensalism relationship 
number_commensalism <-
  sum(method1$Ecological_Relationship == 'Commensalism')


mutualism <-
  method1[method1$Ecological_Relationship == 'Mutualism',]

competition <-
  method1[method1$Ecological_Relationship == 'Competition',]
