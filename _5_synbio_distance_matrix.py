import os as os
import re
import numpy as np
import pandas as pd
import scipy as scipy
import sklearn
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics.pairwise import pairwise_distances
from _6_r_data_processing import to_r_data_processing

#complete
#compare_mat = read_in_binary_matrix(big_matrix_name, location)
def read_in_binary_matrix(ec_binary, name):
    ## Read in the large matrix with the rows and headers as text, all of the data as the binary flags of whether
    ## the EC is present on the genome or not from the Diamond search
    ## Inputs required are the synbio binary summary results and the binary summary matrix for both archaea and bacteria
    ## (*) binary summary for both archaea and bacteria is compile in 0_maintanence script
    ##  Output is the overall summary matrix which includes the combined summary matrix for all (synbio, bacteria, archaea)
    os.chdir(ec_binary)
    binary_matrix_table = pd.read_csv(ec_binary, delimiter=" ", header = 0, index_col = 0)
    print(name, " size of ", np.shape(binary_matrix_table), " successfully imported.")
    ## Opens the matrix that includes bacteria and archaea summary result
    bacteria_binary = pd.read_csv("/home/anna/PycharmProjects/HazID/binary_for_archaea_bacteria.csv",
                                  delimiter="\t", header = 0, index_col = 0)
    synbio_bacteria = bacteria_binary.append(binary_matrix_table)
    #Verically adding synbio matrix to the overall matrix in order to complete comparison, therefore last index
    #should represent the synbio genome that is tested and returned as distances_mat[-1,:] in pass_to_distance.
    #Note, headers are lost and need to be directly passed
    print("Shape of combined matrix using append is: ", np.shape(synbio_bacteria))
    print("Merging of bacteria data and synbio data is complete")
    doc_name_1 = 'complete_binary_matrix.txt'
    synbio_bacteria.to_csv(doc_name_1, header = True, index = True, sep = '\t')
    return synbio_bacteria, doc_name_1

#complete
#    ec_pref = ec_weighting_preference(location)
def ec_weighting_preference(desired_loc):
##TURN ON AFTER COMPLETED TESTING
## Weighting preference should not be used in EC frequency analysis
## Input required is the desired location for the file to save
## Outputs include preference on weighting (Y/N) and the name of the sheet used for scaling
    #print("Would you like to use/upload EC filter/weight scores: (Y/N) ")
    #preference = input()
    preference = "N"
    if preference=="Y":
        print("Would you like to upload your own filter/weight scores? ")
        personal = input()
        if personal == "Y":
            print("Enter the name of your document: " )
            filter_used = input()
            personal_filter_path = os.path.abspath(filter_used)
            if personal_filter_path != desired_loc:
                os.shutil(filter_used, desired_loc)
                default = 0
        elif personal == "N":
            print("You have selected to use the default filter: /home/anna/PycharmProjects/HazID"
                  "/2022_5_24_Proposed_Filter.csv")
            filter_used = "/home/anna/PycharmProjects/HazID/2022_5_24_Proposed_Filter.csv"
            default  = 0
        else:
            print("Error, try again!")
            filter_used = "ERROR"
    else:
        print("No EC filter/weight scores has been selected. All EC numbers will have equal weight")
        default = 1
        filter_used = " "
    #default = 0
    #filter_used = "2022_5_24_Proposed_Filter.csv"
    return default, filter_used

#complete
## Takes in the overall summary matrix, user-inputted settings, as well as the name of document used for EC weights
## Returns the EC space matrix with or without scaling, and the name of genomes found as index names
#    ec_results = ec_weight_implementation(compare_mat, setting, filter_name)
def ec_weight_implementation(synbio_bacteria, pref_set, ec_name):
    #saves the names of genomes in a data frame
    genome_names = pd.DataFrame(synbio_bacteria.index)
    print(genome_names)
    # implements EC filter on the big dataset. Multiplies all of the columns by the EC weight
    print("Beginning EC analysis")
    if pref_set == 0:
        ec_filter = pd.read_csv(ec_name, delimiter = "\t", header = 0, index_col = 0)
        #Transposes EC matrix so the EC scoring will be horizontal
        ec_transposed = ec_filter.T
        #Multiplies values without reading the headers
        input_df = pd.np.multiply(synbio_bacteria, ec_transposed)
    else:
        print("Overall summary matrix remains unchanged")
        input_df = synbio_bacteria
        #print(input_df)
    #synbio_bacteria= pd.concat([genome_names, input_df])
    #np.savetxt("weighted_ec_space.txt", input_df)
    return input_df, genome_names

#takes in the taxonomic names processed in 0_maintenance script
## Bins the raw data into user-specified rank, and averages the EC scores per binning
## Inputs require the EC space as well as the synbio ID for file creation
## Outputs are the EC spaced clustered that will be used for the distance calculation
#   clustered = tax_clustering(ec_space_editted, genome_ID)
def tax_clustering(ec_space, synbio):
    #All of the column names that include '.' are selected (only the EC 1 and 0)
    ec_col = [col for col in ec_space.columns if "." in col]
    full_lineage =  pd.read_csv("/home/anna/PycharmProjects/HazID/taxonomy.tsv", delimiter  = '\t', header = 0)
        #number of rows are not even, code will not work otherwise
    #Change to user inputting taxid
    # print("Enter the known taxonomy of your synbio organism. If no data is available, list SynBio as an entry: ")
    # print("Kingdom: ") ; kingdom = input()
    # print("Phylum: ") ; phylum = input()
    # print("Class: ") ; clas = input()
    # print("Order: ") ; order = input()
    # print("Family: ") ; family = input()
    # print("Genus: ") ; genus = input()
    ##TURN OFF OF TESTING
    kingdom = "Synbio"
    phylum = "Synbio"
    clas= "Synbio"
    order = "Synbio"
    family = "Synbio"
    genus = "Synbio"
    #Creates a dataframe with the lineage of the synbio unit to append to the taxonomy document
    blank_entries = { "Kingdom" : kingdom,
                      "Phylum" : phylum,
                      "Class": clas,
                      "Order": order,
                      "Family": family,
                      "Genus": genus
                      }
    # Appending synbio dataframe to the bottom of the full_lineage dataframe
    full_lineage = full_lineage.append(blank_entries, ignore_index = True)
    print("Lineage matrix size for the complete genomes is: ", full_lineage.shape)
    # concating full lineage df with the binary bacterial one, axis=1 means side by side
    #ec_space_lin = pd.concat([full_lineage.reset_index(drop=True), ec_space.reset_index(drop=True)], axis= 1)
    #print("Would you like to cluster based on lineage? Y/N: ")
    #lineage_preference = input()
    lineage_preference = "Y"
    if lineage_preference == "Y":
        print("Enter the rank you would like to sort by: (Kingdom, Phylum, Class, Order, Family, Genus)")
        cluster_rank = input()
        #This isolates the column that will be used for rank binning
        corresponding_rank = pd.DataFrame(full_lineage[cluster_rank], dtype = "string")
        corresponding_rank.columns = [cluster_rank]
        ec_space_lin = pd.concat([corresponding_rank.reset_index(drop=True),
                                  ec_space.reset_index(drop=True)], axis= 1)
        #np.savetxt('phylum_and_ecspace.txt', ec_space_lin, fmt='%s')
        # grouping the same rows based on a common column, all of the ec scores will be averaged together based on their
        # column value
        ec_grouped = ec_space_lin.groupby(cluster_rank, as_index=True).mean()
        #sets the selected rank as index
        #ec_grouped.set_index(cluster_rank, inplace = True, drop = True)
        #creates the output data frame for distance calculations
        index_names_of_cluster  = pd.DataFrame(ec_grouped.index)
    else:
        #No clustering options means that the dataset will not be binned based on rank and all samples will be calculated
        ec_grouped = ec_space
        cluster_rank = "noClustering"
        index_names_of_cluster = pd.DataFrame()
    ec_grouped = pd.DataFrame(ec_grouped)
    doc_name_2 = "ec_space_" + cluster_rank.lower() + ".txt"
    ec_grouped.to_csv(doc_name_2, header=True, index=True, sep='\t')
    #np.savetxt("ec_space_"+cluster_rank.lower()+".txt", ec_grouped, delimiter = ' ')
    return ec_grouped, lineage_preference, index_names_of_cluster, cluster_rank, doc_name_2

#complete
## Inputs required are the EC space data frame that has been weighted/clustered based on user preference, genome names
## in case no option has been selected, the clustering preference used to consdense data (Y/N), and the index names if
# clustering was chosen
#    synbio_clustered_distances = calculating_distance(clustered_ec, list_genomes, clustering_pref, index_names, genome_ID)
def calculating_distance(input_df, genome_names, lineage_preference, index_names_of_cluster, genome_ID, rank):
    ## Calculate the distances of the datafile using the pdist funciton from scipy. The intial return will only
    ## provide the upper half of the comparisons (row-wise), to crease a symetrical matrix, then create squareform
    doc_name_3 = 'Distance_Matrix_Combined.txt'
    distances_parallel = pairwise_distances(X = input_df, metric = 'euclidean', n_jobs = 18)
    print("Shape of the entire distance matrix is: ", np.shape(distances_parallel))
    ## NaNs at this stage may indicate header name mismatch - check if EC numbers are aligning
    ## If wanting to save the full distance matrix, turn the following flag on. Note: Well above 10 GB
    distances_parallel = pd.DataFrame(distances_parallel)
    distances_parallel.to_csv(doc_name_3, header=False, index= False, sep = '\t')
    if lineage_preference=="Y":
        #If the user chose to cluster, distance matrix will have the rank-related index names for synbio
        distances_synbio = pd.concat([index_names_of_cluster.reset_index(drop=True),
                                      distances_parallel.reset_index(drop=True)], axis= 1)
        #Resets index so the rank column will be the top
        distances_synbio.set_index(rank, inplace=True, drop=True)
        #Finds the row that contains the Synbio unit
        synbio_row = distances_synbio.loc['Synbio', :]
        synbio_column = synbio_row.T
        # Creates the horizontal genome names for the resulting matrix
        synbio_column.index = index_names_of_cluster
    else:
        # Converts the distance matrix of synbio as a data frame. Concat binds the row names data frame with the synbio distance
        # matrix. Result should be a [2,#of total genomes]
        #df_parallel = pd.DataFrame(distances_parallel).reset_index(drop=True)
        distances_synbio = pd.concat([genome_names.reset_index(drop=True),
                                      distances_parallel.reset_index(drop=True)], axis= 1)
        distances_synbio.set_index('Name_of_Genome', inplace = True, drop = True)
        synbio_row = distances_synbio.loc[genome_ID, :]
        synbio_column = synbio_row.T
        synbio_column.index = genome_names
    synbio_column.columns = ['Genome_Name', 'Calculated Distance']
    return synbio_column, doc_name_3

# distance_list_for_synbio = pass_to_distance(name, sb_name, desired_location)
def pass_to_distance(big_matrix_loc, genome_ID, location):
    ## Grabbing the big matrix name and the genome ID that the user provides.
    ## Returns a list of distances found and the name of the document used to save
    [compare_mat, doc_named_1] = read_in_binary_matrix(big_matrix_loc, genome_ID)
    ## asks if you want to apply an EC filter, default is to have all ECs weighted equally. Returns filter preference
    ##  and the name of filter document
    ec_pref = ec_weighting_preference(location)
    ## default or applied filter
    setting = ec_pref[0]
    ## name of the filter document
    filter_name = ec_pref[1]
    ## Transposes EC filter in order to have same dimensions as the binary EC space
    # Completes column wise matrix multiplication, returns the weighted EC space
    ec_results = ec_weight_implementation(compare_mat, setting, filter_name)
    ec_space_editted = ec_results[0]
    list_genomes = ec_results[1]
    ## Takes in the weighted EC space and breaks down the name of genome column into genus and species takes average of
    # ec values and combines genomes into a genus based data frame returns the clustered matrix of ec scores
    clustered = tax_clustering(ec_space_editted, genome_ID)
    clustered_ec = clustered[0]
    clustering_pref = clustered[1]
    index_names = clustered[2]
    rank = clustered[3]
    doc_name_2 = clustered[4]
    ## Completes distance calculation on an euclidean basis. Returns the synbio column with genome names appended
    [synbio_clustered_distances, doc_name_3] = calculating_distance(clustered_ec, list_genomes, clustering_pref, index_names, genome_ID, rank)
    #Saves the output from distance matrix
    name= genome_ID+"_"+rank+"_ClusteredDistanceMatrix.txt"
    #Saves only the synbio column
    synbio_clustered_distances.to_csv(name, header=True, index=True, sep = '\t')
    # Sends documents for R data processing...
    to_r_data_processing(doc_named_1,doc_name_2, doc_name_3, location)
    return synbio_clustered_distances

