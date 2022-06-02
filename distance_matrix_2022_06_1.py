import os as os
import re
import numpy as np
import pandas as pd
import scipy as scipy
import sklearn
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics.pairwise import pairwise_distances

# This function completes a genome to genome comparative analysis.

#complete
def read_in_binary_matrix(name):
    ## Read in the large matrix with the rows and headers as text, all of the data as the binary flags of whether
    ## the EC is present on the genome or not from the Diamond search
    binary_matrix_table = pd.read_csv(name, delimiter=" ", header = 0, index_col = 0)
    print(name, " size of ", np.shape(binary_matrix_table), " successfully imported.")
    bacteria_binary = pd.read_csv("bacteria_big_matrix_tab.txt", delimiter="\t", header = 0, index_col = 0)
    synbio_bacteria = bacteria_binary.append(binary_matrix_table)
    #Verically adding synbio matrix to the overall matrix in order to complete comparison, therefore last index
    #should represent the synbio genome that is tested and returned as distances_mat[-1,:] in pass_to_distance.
    #Note, headers are lost and need to be directly passed
    print("Shape of combined matrix using append is: ", np.shape(synbio_bacteria))
    print("Merging of bacteria data and synbio data is complete")
    print(type(synbio_bacteria))

    return synbio_bacteria

#complete
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
        input_df =  synbio_bacteria
    #synbio_bacteria= pd.concat([genome_names, input_df])
    #np.savetxt("weighted_ec_space.txt", input_df)
    return input_df, genome_names

#takes in the taxonomic names processed in R
def tax_clustering(ec_space, synbio):
    ec_col = [col for col in ec_space.columns if "." in col]  # all columns that contain the ones and zeros
    full_lineage =  pd.read_csv("taxonomy_clustering.txt", delimiter  = '\t', header = 0)        #number of rows are not even, code will not work otherwise
    #Change to user inputting taxid
    #print("Enter the known taxonomy of your synbio organism: ")
    #print("Kingdom: ") kingdom = input()
    #print("Phylum: ") phylum = input()
    # print("Class: ") clas = input()
    # print("Order: ") order = input()
    # print("Family: ") family = input()
    # print("Genus: ") genus = input()
    kingdom = "Bacteria"
    phylum = "Proteobacteria"
    clas= "Gammaproteobacteria"
    order = "Enterobacterales"
    family = "Enterobacteriaceae"
    genus = "Escherichia"
    blank_entries = { "Kingdom" : kingdom,
                      "Phylum" : phylum,
                      "Class": clas,
                      "Order": order,
                      "Family": family,
                      "Genus": genus
                      }
    full_lineage = full_lineage.append(blank_entries, ignore_index=True)
    print("Lineage matrix size for the complete genomes is: ",full_lineage.shape)
    ec_space_lin = pd.concat([full_lineage.reset_index(drop=True), ec_space.reset_index(drop=True)], axis= 1)         #concating full lineage df with the binary bacterial one, axis=1 means side by side
    print("Would you like to cluster based on lineage? Y/N: ")
    lineage_preference = input()
    if lineage_preference == "Y":
        print("Enter the rank you would like to sort by: (Kingdom, Phylum, Class, Order, Family, Genus)")
        cluster_rank = input()
        ec_grouped = ec_space_lin.groupby(cluster_rank, as_index= False)[ec_col].agg('mean')         #grouping the same rows based on a common column
        #sets the phylum as index
        ec_grouped.set_index(cluster_rank, inplace = True, drop = True)
        index_names_of_cluster  = pd.DataFrame(ec_grouped.index)
    else:
        ec_grouped = ec_space
        #ec_grouped.set_index("Name_of_Genome", inplace = True, drop= True)
        index_names_of_cluster = pd.DataFrame()

    np.savetxt("ec_space_"+cluster_rank.lower()+".txt", ec_grouped)
    return ec_grouped, lineage_preference, index_names_of_cluster

#complete
def calculating_distance(input_df, genome_names, lineage_preference, index_names_of_cluster):
    ## Calculate the distances of the datafile using the pdist funciton from scipy. The intial return will only
    ## provide the upper half of the comparisons (row-wise), to crease a symetrical matrix, then create squareform
    distances_parallel = pairwise_distances(X = input_df, metric = 'euclidean', n_jobs = 18)
    ## NaNs at this stage may indicate header name mismatch - check if EC numbers are aligning
    ## If wanting to save the full distance matrix, turn the following flag on. Note: Well above 10 GBs
    #np.savetxt("Distance_matrix_Combined.txt",distances_parallel)
    #print(type(distances_parallel))
    if lineage_preference=="Y":
         distances_synbio = pd.concat([index_names_of_cluster.reset_index(drop=True), pd.DataFrame(distances_parallel).reset_index(drop=True)], axis= 1)
    else:
        # Converts the distance matrix of synbio as a data frame. Concat binds the row names data frame with the synbio distance
        # matrix. Result should be a [2,#of total genomes]
        #df_parallel = pd.DataFrame(distances_parallel).reset_index(drop=True)
        distances_synbio = pd.concat([genome_names.reset_index(drop=True), pd.DataFrame(distances_parallel).reset_index(drop=True)], axis= 1)
        #distances_synbio.columns = ['Genome_Name', 'Calculated Distance']
    return distances_synbio

#complete
def ec_weighting_preference(desired_loc):
##TURN ON AFTER COMPLETED TESTING
    # print("Would you like to use/upload EC filter/weight scores: (Y/N) ")
    # preference = input()
    # if preference=="Y":
    #     print("Would you like to upload your own filter/weight scores? ")
    #     personal = input()
    #     if personal == "Y":
    #         print("Enter the name of your document: " )
    #         filter_used = input()
    #         personal_filter_path = os.path.abspath(filter_used)
    #         if personal_filter_path != desired_loc:
    #             os.shutil(filter_used, desired_loc)
    #             default = 0
    #     elif personal == "N":
    #         print("You have selected to use the default filter: 2022_5_24_Proposed_Filter.csv")
    #         filter_used = "2022_5_24_Proposed_Filter.csv"
    #         default  = 0
    #     else:
    #         print("Error, try again!")
    #         filter_used = "ERROR"
    # else:
    #     print("No EC filter/weight scores has been selected. All EC numbers will have equal weight")
    #     default = 1
    #     filter_used = " "
    default = 0
    filter_used = "2022_5_24_Proposed_Filter.csv"
    return default, filter_used


def pass_to_distance(big_matrix_name, genome_ID, location):
    ## Grabbing the big matrix name and the genome ID that the user provides.
    ## Returns a list of distances found
    compare_mat = read_in_binary_matrix(big_matrix_name)
    ## asks if you want to apply an EC filter, default is to have all ECs weighted equally. Returns filter preference
    ##  and the name of filter document
    ec_pref = ec_weighting_preference(location)
    ## default or applied filter
    setting = ec_pref[0]
    ## name of the filter document
    filter_name = ec_pref[1]
    ## Transposes EC filter in order to have same dimensions as the binary EC space
    ## Completes column wise matrix multiplication, returns the weighted EC space
    ec_results = ec_weight_implementation(compare_mat, setting, filter_name)
    ec_space_editted = ec_results[0]
    list_genomes = ec_results[1]
    ## Takes in the weighted EC space and breaks down the name of genome column into genus and species takes average of ec values and combines genomes into a genus based data frame
    ## returns the clustered matrix of ec scores
    clustered = tax_clustering(ec_space_editted, genome_ID)
    clustered_ec = clustered[0]
    clustering_pref = clustered[1]
    index_names = clustered[2]
    ## Completes distance calculation on an euclidean basis. Returns the synbio column with genome names appended
    synbio_clustered_distances = calculating_distance(clustered_ec, list_genomes, clustering_pref, index_names)
    return synbio_clustered_distances


