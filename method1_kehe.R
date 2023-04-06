#install.packages(c('stringr'))
library(stringr)
library(janitor)
library(data.table)
# method 1
source('scoring.R')
################################################################################

## Overview: Finds 1000 random pairs to test, constructs a binary EC overlap matrix, based on whether
## the EC number is present in both pairs (=1), finds EC overlap frequency
## Requires the EC BSM for each genome in a cluster

m1_ec_overlap_btwn_sets <- function(bug1, bug2, b1, b2) {
  # Creates a dataframe with 1000 random samples for Cluster 1
  bug1_replicate_sample <-
    data.frame(replicate(1000, sample(bug1$Name_of_Genome)))
  # Isolates the row with the random genomes selected for Cluster 1
  bug1_random_1000 <- transpose(bug1_replicate_sample[2, ])
  # Creates a dataframe with 1000 random samples for Cluster 2
  bug2_replicate_sample <-
    data.frame(replicate(1000, sample(bug2$Name_of_Genome)))
  # Isolates the row with the random genomes selected for Cluster 2
  bug2_random_1000 <- transpose(bug2_replicate_sample[2, ])
  colnames(bug1_random_1000) <- c('Name_of_Genome')
  colnames(bug2_random_1000) <- c('Name_of_Genome')
  # Merges randomized genome names data frame with the EC BSM, finding the EC BSM
  # for 1000 randomly picked genomes for Cluster 1
  b1_random_bsm <- merge(bug1_random_1000, bug1, by = 'Name_of_Genome')
  # Merges randomized genome names data frame with the EC BSM, finding the EC BSM
  # for 1000 randomly picked genomes for Cluster 2
  b2_random_bsm <- merge(bug2_random_1000, bug2, by = 'Name_of_Genome')
  
  ## Overview: Randomly picked genomes from each cluster are compared to each other. If both have the same EC number,
  ## a new dataframe (combos_added) will have a value of 1, If one genome has the EC number but the other one doesn't, the
  ## new dataframe will have a value of zero
  
  # Creates row names for cluster-to-cluster comparison for combos_added dataframe which will contain the EC comparisons
  combined_row_names <-
    data.frame(paste(
      b1_random_bsm$Name_of_Genome,
      b2_random_bsm$Name_of_Genome,
      sep = "|"
    ))
  colnames(combined_row_names) <- "Genome_Pairs"
  # For any repeated pairs, the names is modified to contain the replicate number
  combined_row_names$Genome_Pairs <-
    make.unique(combined_row_names$Genome_Pairs, sep = '_replicate_')
  combined_row_names <- data.frame(combined_row_names)
  # Isolates the binary field for Cluster 1
  bug1_binary <- b1_random_bsm[, 3:8198]
  # Isolates the binary field for Cluster 2
  bug2_binary <- b2_random_bsm[, 3:8198]
  # Completes addition between Cluster 1 and Cluster 2 fields
  combos_added <- bug1_binary + bug2_binary
  # Changes 1 to a zero
  combos_added[combos_added == 1] <- 0
  # Changes 2 to a 1
  combos_added[combos_added == 2] <- 1
  # Sets the first row as column names and removes the first row
  names(combined_row_names) <- NULL
  
  row.names(combos_added) <- combined_row_names[, 1]
  # Adds the EC frequencies together and divides by the total genomes sampled (=1000)
  # Creates a matrix of EC numbers vs. average (1 row * number of ECs)
  ec_overlap_average <-
    data.frame(colSums(combos_added) / nrow(combos_added))
  # Transposes matrix so it becomes a column of EC frequencies
  ec_overlap_frequency <- transpose(ec_overlap_average)
  row.names(ec_overlap_frequency) <- ('EC_Overlap_Frequency')
  # Adds EC numbers as row names
  colnames(ec_overlap_frequency) <- colnames(combos_added)
  # Adds labels to a dataframe which has the EC Overlap frequencies and the binary field
  binary_w_frequency <- rbind(combos_added, ec_overlap_frequency)
  #write.csv(binary_w_frequency,'/home/anna/rstudio-2022.12.0-353-amd64/Friedman/ec_overlpap_frequencies.txt')
  to_next_function <-
    data.frame(rbind(
      colnames(binary_w_frequency),
      tail(binary_w_frequency, n = 1)
    ))
  to_next_function <- data.frame(transpose(to_next_function))
  colnames(to_next_function) <- c('EC', 'Frequency')
  write.csv(
    to_next_function,
    '/home/anna/rstudio-2022.12.0-353-amd64/Friedman/ec_overlap_sent_to_next_function.txt'
  )
  inchi_key_lineup(to_next_function, b1, b2)
}
################################################################################

## Overview: Creates a binary presence matrix of InChI-Keys per EC Number
## Adds the InChI-Key presences per EC and multiplies by the frequency
## Score is assigned by summing all of the product found above

# Competition
inchi_key_lineup <- function(ec_overlap, b1, b2) {
  ec_overlap$EC <- gsub('X', '', ec_overlap$EC)
  ec_numbers <- ec_overlap$EC
  ec_frequency <- ec_overlap[, 2]
  # List of all InChI-Keys publicly available was downloaded from MetaCyc Database on 2/9/2023
  list_of_all_inchikeys <-
    read.csv(
      '/home/anna/rstudio-2022.12.0-353-amd64/Friedman/all_inchikeys.txt',
      header = TRUE,
      quote = '',
      sep = '\t'
    )
  colnames(list_of_all_inchikeys) <- c('Compounds', 'InChI_Key')
  # Opens a document with all reactions and all InChI-Keys associated
  
  metacyc_inchi_key_list <-
    read.csv(
      '/home/anna/Desktop/EcoGenoRisk/HazID/CompetitorFind/All-reactions-of-MetaCyc.txt',
      header = TRUE,
      quote = '',
      sep = '\t'
    )
  colnames(metacyc_inchi_key_list) <-
    c(
      "Reaction",
      "EC",
      "Substrates",
      "Substrate_InChI_Key",
      "Reactants_of_reaction",
      "Reactants_InChI_Key",
      "Products_of_reaction",
      "Products-InChI_Key"
    )
  inchi_key <- data.frame(metacyc_inchi_key_list$Reactants_InChI_Key)
  colnames(inchi_key) <- c('Reactants_InChI_Key')
  ec <- data.frame(metacyc_inchi_key_list$EC)
  colnames(ec) <- c('EC')
  # Subsets the MetaCyc dataframe to contain only the EC numbers and InChI-Key list
  ec_inchikey <- data.frame(ec, inchi_key)
  colnames(ec_inchikey) <- c('EC', 'Reactants_InChI_Key')
  # Turn on to save MetaCyc dataframe subset
  #write.csv(ec_inchikey,'/home/anna/rstudio-2022.12.0-353-amd64/Friedman/ec_inchikey_subdf.txt')
  
  # Creates dimensions for an empty data frame
  columns <- nrow(list_of_all_inchikeys)
  rows <- nrow(ec_overlap)
  # Creates an empty data frame which will contain the InChI-Key presence per EC number
  blank <- data.frame(matrix(
    data = 0,
    ncol = columns,
    nrow = rows
  ))
  # Creates the the InChI-Key presence matrix with labelled EC numbers, EC overlap frequency along the vertical axis
  inchikey_present <- cbind(ec_numbers, ec_overlap, blank)
  # Isolates a data frame of InChI-Keys only, adds a column name
  list_of_all_inchikey <- data.frame(list_of_all_inchikeys$InChI_Key)
  colnames(list_of_all_inchikey) <- c('InChI_Key')
  # Removes any blank rows from the InChI-key list
  list_of_all_inchikey <-
    data.frame(list_of_all_inchikey[-which(list_of_all_inchikey$InChI_Key ==
                                             ""), ])
  colnames(list_of_all_inchikey) <- c('InChI_Key')
  # Adds column names to the blank data frame
  colnames(inchikey_present) <-
    cbind('EC', 'EC_Overlap_Frequency', t(list_of_all_inchikey))
  # Removes unnecessary columns
  inchikey_present <- inchikey_present[, -1]
  row.names(inchikey_present) <- ec_numbers
  
  ## Overview: Cycles through EC list one by one
  ## Using the MetaCyc subset dataframe, isolates those rows that include the EC number
  ## Cycles through the isolated row, and sees if there is an InChI-Key match 
  ## If there is an InChI-Key match, enters a 1 in the appropriate place in the dataframe
  ## Note: extremely slow, needs reworking 
  
  # Selects an EC from an EC number list
  # Initialize the inchikey_present matrix
  inchikey_present <- matrix(0, nrow = nrow(ec_overlap), ncol = nrow(list_of_all_inchikey))
  
  # Loop over rows of the ec_inchikey data frame
  apply(ec_inchikey, 1, function(x) {
    # Find the rows in the MetaCyc subset dataframe where the EC number is found
    ec_rows <- ec_inchikey[grep(x[1], ec_inchikey$EC), ]
    # Check if there are any recorded data for the EC number
    if (nrow(ec_rows) > 0) {
      # Loop over the InChI_Keys
      for (j in 1:nrow(list_of_all_inchikey)) {
        # Check if the InChI_Key is present in the Reactants_InChI_Key column
        if (any(str_detect(x[2], list_of_all_inchikey$InChI_Key[j]))) {
          # Mark the corresponding element in the inchikey_present matrix as 1
          inchikey_present[which(ec_overlap$EC == ec_rows$EC), j] <- 1
        }
      }
    }
  })
  
  # Turn on to see the dimensions of the blank dataframe
  # print(dim(inchikey_present))
  # Saves the blank dataframe as a text file
  write.csv(
    inchikey_present,
    '/home/anna/rstudio-2022.12.0-353-amd64/Friedman/inchikey_present.txt'
  )
  # Adds the InChI-Key presence in the data frame and multiplies by the EC overlap
  ec_frequency_inchikey_presence <-
    rowSums(list_of_all_inchikeys[, 1:18788]) * list_of_all_inchikeys$EC_Overlap_Frequency
  # Adds the product together to assign a final score
  final_score <- colSums(ec_frequency_inchikey_presence)
  # Sends the final score into the scoring function on the scoring.R file for formatting
  scoring(final_score, b1, b2, 'method1')
}
