#This is the _4_ and _5_ linearized, for conceptual understanding and to streamline the code when the time comes
##--------------------------------------------------------------------------------------------#
#[distance_list_for_synbio, new_loc ]= pass_to_distance(synbio_binary, sb_name, desired_location)
#sb_name = 'GCF_002926195.1_ASM292619v1_protein'
#synbio_binary = analysis_output[0]
#desired_location = desired_location + "/" + sb_name

import pandas as pd 
import numpy as np #In this case, /projects/jodo9280/EcoDr/GCF_002926195.1_ASM292619v1_protein/DIAMOND_matches/GCF_002926195.1_ASM292619v1_protein
import os
from sklearn.metrics import pairwise_distances
##------------------------------------------------------------------------------------------------------------------------------#
def pass_to_distance(synbio_binary, sb_name, desired_location):
    [compare_mat, doc_named_1] = read_in_binary_matrix(synbio_binary, sb_name)
    ec_pref = ec_weighting_preference(desired_location)
    # Receives preference for default or applied filter
    setting = ec_pref[0]
    # Returns name of the filter document
    filter_name = ec_pref[1]
    # Transposes EC filter to have same dimensions as the binary EC space
    # Completes column-wise matrix multiplication, returns the weighted EC space
    ec_results = ec_weight_implementation(compare_mat, setting, filter_name)
    ec_space_editted = ec_results[0]
    print('EC_spec_editted')
    print(ec_space_editted)
    list_genomes = ec_results[1]
    # Takes in the weighted EC space and breaks down the name of genome column into genus and species takes average of
    # EC values and combines genomes into a genus based data frame returns the clustered matrix of EC scores
    clustered = tax_clustering(ec_space_editted, sb_name)
    clustered_ec = clustered[0]
    clustering_pref = clustered[1]
    index_names = clustered[2]
    rank = clustered[3]
    doc_name_2 = clustered[4]
    doc_name_3 = clustered[5]
    # Completes distance calculation on an Euclidean basis. Returns the synbio column with genome names appended
    [synbio_clustered_distances, name] = calculating_distance(clustered_ec, list_genomes, clustering_pref,
                                                              index_names, sb_name, rank)
    # Saves the output from distance matrix
    #doc_name_4 = name
    # Saves only the synbio column
    # Saves only the synbio column
    synbio_clustered_distances.to_csv(name, header=True, index=True, sep='\t')
    return synbio_clustered_distances, desired_location
##------------------------------------------------------------------------------------------------------------------------------#
def read_in_binary_matrix(synbio_binary, sb_name):
    # Converts synbio summary matrix into a dataframe
    synbio_binary = pd.read_csv(synbio_binary, delimiter=" ", header=0)
    print(synbio_binary)
    synbio_binary = synbio_binary.set_index('Name_of_Genome')
    #index is a label for all rows - allows the two seperate dataframes to come together since they share an index
    print(synbio_binary)
    print(sb_name, " size of ", np.shape(synbio_binary), " successfully imported.")
    # Opens the matrix that includes the Bacteria and Archaea summary result
    domain_binary = pd.read_csv('/projects/jodo9280/EcoDr/EcoDr/EcoDr_binary_matrix.txt',
                                  delimiter=" ", header=0)
    #this is from chunk 1, and is the EC_Binary we generated earlier
    domain_binary = domain_binary.set_index('Name_of_Genome')
    print(domain_binary)
    # Sends to a function for direct genome to genome comparison based on EC summary matrix
    genome_to_genome_diffcomp(synbio_binary, domain_binary)
    synbio_bacteria = pd.concat([domain_binary, synbio_binary])
    # Verically adding synbio matrix to the overall matrix to complete the comparison, therefore last index
    # should represent the synbio genome that is tested and returned as distances_mat[-1,:] in pass_to_distance.
    # Note, headers are lost and need to be directly passed
    print("Shape of combined matrix using append is: ", np.shape(synbio_bacteria))
    print("Merging of bacteria data and synbio data is complete")
    doc_name_1 = 'complete_binary_matrix.txt'
    # Removes any duplicates
    synbio_bacteria = synbio_bacteria[~synbio_bacteria.index.duplicated(keep='first')]
    synbio_bacteria.to_csv(doc_name_1, header=True, index=True, sep='\t')
    return synbio_bacteria, doc_name_1
##-------------------------------------------------------------------------------------------------------------------------##
def genome_to_genome_diffcomp(synbio_binary, domain_binary):
    names_of_orgs = pd.DataFrame(domain_binary.index) 
    diff = pd.DataFrame(abs(domain_binary.values - synbio_binary.values))
    row_sum = diff.sum(axis=1)
    df1 = pd.DataFrame(row_sum)
    difference_based_comparison = pd.concat([names_of_orgs, df1], axis=1)
    difference_based_comparison.columns = ['Organisms Compared to Synbio', 'Difference Score']
    difference_based_comparison = difference_based_comparison.sort_values(by='Difference Score',
                                                                          ignore_index=True).reset_index(drop=True)
    difference_based_comparison.to_csv('Difference_Based_Comparison_Score.txt', header=True, index=True, sep='\t')
##----------------------------------------------------------------------------------------------------------------------------###
#prioritizes certain EC categories
#need to track down the file here 
#weighting will change the binary so instead of a one you may have a 10 
def ec_weighting_preference(desired_loc):
    print("Would you like to use/upload EC filter/weight scores: (Y/N) ")
    # preference = input()
    # Fix this later, turned off 20231117 by CBM
    preference = "N"
    if preference == "Y":
        print("Would you like to upload your own filter/weight scores? ")
        personal = input()
        if personal == "Y":
            print("Enter the filepath of your document: ")
            filter_used = input()
            personal_filter_path = os.path.abspath(filter_used)
            if personal_filter_path != desired_loc:
                os.shutil(filter_used, desired_loc)
                default = 0
        elif personal == "N": #inport in EC list where everything is 1 as a default filter 
            print("You have selected to use the default filter: /home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap"
                  "/2022_5_24_Proposed_Filter.csv")
            filter_used = "/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/2022_5_24_Proposed_Filter.csv"
            default = 0
        else:
            print("Error, try again!")
            filter_used = "ERROR"
    else:
        print("No EC filter/weight scores has been selected. All EC numbers will have equal weight")
        default = 1
        filter_used = " "
    return default, filter_used
##---------------------------------------------------------------------------------------------------------------------------------###
def ec_weight_implementation(synbio_bacteria, pref_set, ec_name):
    # Saves the names of genomes in a data frame
    genome_names = pd.DataFrame(synbio_bacteria.index)
    print(genome_names)
    # Implements EC filter on the big dataset. Multiplies all of the columns by the EC weight
    print("Beginning EC analysis")
    if pref_set == 0:
        ec_filter = pd.read_csv(ec_name, delimiter="\t", header=0, index_col=0)
        # Transposes the EC matrix to create a horizontal EC weighted matrix
        ec_transposed = ec_filter.transpose()
        # Multiplies original ec space by the chosen filter without reading the headers
        input_df = pd.np.multiply(synbio_bacteria, ec_transposed)
        ec_filter = 'weighted'
    else:
        print("Overall summary matrix remains unchanged")
        input_df = synbio_bacteria
        ec_filter = 'unweighted'
    # Turn on to save the weighted EC scores
    np.savetxt(ec_filter +"_ec_space.txt", input_df)
    return input_df, genome_names
##-------------------------------------------------------------------------------------------------------------------------------------###
def tax_clustering(ec_space, synbio):
    # All of the column names that include '.' are selected (only the EC 1 and 0)
    ec_col = [col for col in ec_space.columns if "." in col]
    doc_3_name = '/projects/jodo9280/EcoDr/taxonomy_2023_6_27.tsv' #where does this come from. 
    full_lineage = pd.read_csv(doc_3_name, delimiter='\t',
                               header=0)
    # If the number of rows are not even, then code will not work otherwise
    # Asks if the user wants to use default lineage (synbio entered for all ranks), or it they want to manually enter the taxonomy
    print(
        "Enter the known taxonomy of your synbio organism. If no data is available, list Synbio as an entry for all ranks: ")
    print("Kingdom: ");
    # kingdom = input()
    kingdom = "No"
    print("Phylum: ");
    # phylum = input()
    phylum = "Chance"
    print("Class: ");
    # clas = input()
    clas = "of"
    print("Order: ");
    # order = input()
    order = "this"
    print("Family: ");
    # family = input()
    family = "working"
    print("Genus: ");
    genus = "out"
    # Creates a dataframe with the lineage of the synbio unit to append to the taxonomy document
    blank_entries = {"Kingdom": kingdom,
                     "Phylum": phylum,
                     "Class": clas,
                     "Order": order,
                     "Family": family,
                     "Genus": genus
                     }
    # Appends synbio dataframe to the bottom of the full_lineage dataframe
    print(ec_space.columns)
    print(full_lineage.columns)
    ec_space.to_csv('ec_space_before_merging.txt', header=True, index=True, sep='\t')
    ec_space.index.name = 'Name_of_Genome'
    full_lineage = pd.merge(full_lineage,ec_space, on='Name_of_Genome', how='right')
    print("Lineage matrix size for the complete genomes is: ", full_lineage.shape)
    print("Would you like to cluster based on lineage? Y/N: ")
    lineage_preference = "N" # input()
    if lineage_preference == "Y":
        print("Enter the rank you would like to sort by: (Kingdom, Phylum, Class, Order, Family, Genus)")
        cluster_rank = input()
        # This isolates the column that will be used for rank binning
        corresponding_rank = pd.DataFrame(full_lineage[cluster_rank], dtype="string")
        corresponding_rank.columns = [cluster_rank]
        ec_space_lin = pd.concat([corresponding_rank.reset_index(drop=True),
                                  ec_space.reset_index(drop=True)], axis=1)
        # Concats full lineage df with the binary bacterial one, axis=1 means side by side
        # Turn on if you would like to see the original EC binary summary and the full taxonomy lineage for the genome in question
        # np.savetxt('phylum_and_ecspace.txt', ec_space_lin, fmt='%s')
        # Groups the same rows based on a common column, all of the EC scores will be averaged together based on their
        # column value
        ec_grouped = ec_space_lin.groupby(cluster_rank, as_index=True).mean()
        # Creates the output data frame for distance calculations
        index_names_of_cluster = pd.DataFrame(ec_grouped.index)
    else:
        # No clustering options means that the dataset will not be binned based on rank and all samples will be calculated
        ec_grouped = ec_space
        cluster_rank = "noClustering"
        index_names_of_cluster = pd.DataFrame()
    # Creates a dataframe of the EC averages only
    ec_grouped = pd.DataFrame(ec_grouped)
    # Packages name for the EC averaged space
    doc_name_2 = "ec_grouped_by_" + cluster_rank.lower() + ".txt"
    # Saves the averaged EC space, will be used for PCA analysis in R
    ec_grouped.to_csv(doc_name_2, header=True, index=True, sep='\t')
    return ec_grouped, lineage_preference, index_names_of_cluster, cluster_rank, doc_name_2, doc_3_name
##-------------------------------------------------------------------------------------------------------###
def calculating_distance(input_df, genome_names, lineage_preference, index_names_of_cluster, genome_ID, rank):
    #calculating_distance(clustered_ec, list_genomes, clustering_pref, index_names, sb_name, rank)
    # Calculate the distances of the datafile using the pdist funciton from scipy. The intial return will only
    # provide the upper half of the comparisons (row-wise) to create a symetrical matrix then create squareform
    name = genome_ID + "_" + rank + "_unclustered_" + 'Distance_Matrix_Combined.txt'
    distances_parallel = pairwise_distances(X=input_df, metric='euclidean', n_jobs=18)
    #imported function that takes in the ec, calculates the distance between rows, n_jobs is for large datasets
    print("Shape of the entire distance matrix is: ", np.shape(distances_parallel))
    # NaNs at this stage may indicate header name mismatch - check if EC numbers are aligning
    # If wanting to save the full distance matrix, then turn the following flag on. Note: Well above 10 GB
    distances_parallel = pd.DataFrame(distances_parallel)
    print(distances_parallel)
    #print('Chocolate Sunday')
    # Returns the full distance matrix for dendrogram construction in R script
    #distances_parallel.to_csv('diff_weighted_unclustered_distance_matrix.txt', header=False, index=False, sep='\t')
    #distances_parallel.to_csv('both_chimeras_distance_matrix_clustered.txt', header=False, index=False, sep='\t')
    print("Distance matrix is downloaded")
    if lineage_preference == "Y":
        # If the user chose to cluster, then the distance matrix will have the rank-related index names for synbio
        distances_synbio = pd.concat([index_names_of_cluster.reset_index(drop=True),distances_parallel.reset_index(drop=True)], axis=1)
        # Resets index so the rank column will be the top
        distances_parallel.index = index_names_of_cluster
	#might need to fix this above
        distances_parallel.to_csv('Chimera1_unclustered_unweighted_DM.txt', header = True, index = True, sep='\t')
        # Finds the row that contains the Synbio unit
        tsv_name = genome_ID.replace(".faa", "")
        tsv_name2 = tsv_name + "_matches.tsv"
        print(tsv_name2)
        synbio_row = distances_parallel.loc[tsv_name2 , :]
        synbio_column = synbio_row.T
        # Creates the vertical genome names for the resulting matrix
        synbio_column.index = index_names_of_cluster
    else:
        print(genome_names)
        # Converts the distance matrix of synbio as a data frame. Concat binds the row names dataframe with the synbio distance
        # matrix. Result should be a [2,#of total genomes]
        distances_synbio = pd.concat([genome_names.reset_index(drop=True),
                                       distances_parallel.reset_index(drop=True)], axis=1)
        #axis=1 is the columns. Just adds genome names to the final score output
        distances_parallel.index = genome_names['Name_of_Genome'].tolist()
        distances_parallel.to_csv('Chimera1_unclustered_unweighted_DM.txt', header = True, index = True, sep='\t')
        #distances_parallel.set_index('Name_of_Genome', inplace=True, drop=True)
        # Finds the row that contains the synbio genome based on genome ID
        tsv_name = genome_ID.replace(".faa", "")
        tsv_name2 = tsv_name + "_matches.tsv"
        print(tsv_name2)
        synbio_row = distances_parallel.loc[tsv_name2 ,:]
        #grabs the name of the first GCF using the index
        synbio_column = synbio_row.T
        # Creates the vertical genome names for the resulting matrix
        synbio_column.index = genome_names
    # Sets column names
    synbio_column.columns = ['Genome_Name', 'Calculated Distance']
    return synbio_column, name
