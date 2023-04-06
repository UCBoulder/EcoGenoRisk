setwd("/home/anna/rstudio-2022.12.0-353-amd64/Friedman")
# for package help: https://stackoverflow.com/questions/63390194/package-xxx-was-installed-before-r-4-0-0-please-re-install-it
library(stringr)
library(janitor)
library(data.table)
source('method1_kehe.R')
source("method2_kehe.R")

# Iterates through the friedman pairs and isolates the EC BSM for all genomes
# pertaining to the bug pairs
# Sends to method 1 analysis (m1_ec_overlap_btwn_sets) and method 2 analysis 
# (m2_inchikey_list_per_genome)
pair_designation <- function(friedman_pairs, total_genomes) {
  for (line in 1:nrow(friedman_pairs)) {
    bug1 <- friedman_pairs[line, ][1]
    bug2 <- friedman_pairs[line, ][2]
    
    bug1_ec_bsm <- total_genomes[grep(bug1, total_genomes$Code_name), ]
    bug2_ec_bsm <- total_genomes[grep(bug2, total_genomes$Code_name), ]
    
    if (nrow(bug1_ec_bsm) > 0 & nrow(bug2_ec_bsm) > 0) {
      write.csv(bug1_ec_bsm,'/home/anna/rstudio-2022.12.0-353-amd64/Friedman/bug1_ec_bsm.txt')
      m1_ec_overlap_btwn_sets(bug1_ec_bsm, bug2_ec_bsm,bug1,bug2)
      m2_inchikey_list_per_genome(bug1_ec_bsm,bug2_ec_bsm)
    }
  }
}

################################################################################
# Uploading files: complete EC BSM (for Archaea and Bacteria),
# Friedman supplementary document which lists all of the genome clusters and found relationship
binary_matrix <-
  read.delim(
    'combined_binary_summary_matrix.csv',
    header = TRUE,
    quote = '',
    sep = '\t'
  )
bug_by_genome <-
  read.delim('genome_species_assigned_by_codename.txt',
             header = TRUE,
             quote = '',
             sep = '\t')
friedman_testcases1 <-
  read.delim(
    'kehe_output_summary.csv',
    header = TRUE,
    quote = '',
    sep = ','
  )

# Reformats dataframes 
binary_matrix <- as.data.frame(binary_matrix)
header <- names(binary_matrix)
complete_header <- gsub('X', '', header)
names(binary_matrix) <- complete_header

# Typecasts the binary field as numeric data type
binary_matrix[, 2:8197] <- sapply(binary_matrix[, 2:8197], as.numeric)

# Isolates the genome names per cluster-type (referred to as "bug")
bug_by_genome <-
  as.data.frame(bug_by_genome[, !names(bug_by_genome) %in% c('X', 'Species')])


friedman_testcases <- as.data.frame(friedman_testcases1)
header2 <- names(friedman_testcases)
complete_header2 <- gsub('X', '', header2)
header3 <- str_replace(complete_header2, '.', '')
names(friedman_testcases) <- header3


# Setting up cluster groups
# Finds the EC BSM for each genome identified in a cluster
total_genomes <-
  merge(bug_by_genome, binary_matrix, by = 'Name_of_Genome')
# Removes any empty rows where there is no species specified
total_genomes_noempty <-
  total_genomes[!(total_genomes$Name_of_Genome == ''), ]
# Gets rid of unneccessary columns that were vestiges of previous dataframes
total_genomes_noempty <-
  total_genomes_noempty[, !names(total_genomes_noempty) %in% c('X', 'Species')]

# Selects relationships where Bacteria were grown on sucrose
friedman_glucose <-
  friedman_testcases[friedman_testcases[, 3] == '"Glucose"', ]

# Subsets the friedman dataframe to include only the Bug1 and Bug2 clusters grown on glucose
friedman_glucose_bugs <-
  as.data.frame(cbind(friedman_glucose$Bug.1., friedman_glucose$Bug.2.))
colnames(friedman_glucose_bugs) <- c('Bug1', 'Bug2')
# String reformatting
friedman_glucose_bugs[] <-
  lapply(friedman_glucose_bugs,
         gsub,
         pattern = '"',
         replacement = '')

# Sends to pair designation, this function requires the complete list of cluster pairs
# and the non-empty EC BSM for genomes found in the clusters
pair_designation(friedman_glucose_bugs, total_genomes_noempty)
