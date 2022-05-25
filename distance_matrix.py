import os as os
import re
import numpy as np
import pandas as pd
import scipy as scipy
import sklearn
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics.pairwise import pairwise_distances

# This function completes a genome to genome comparative analysis.


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


def calculate_dsitance(synbio_bacteria, pref_set, ec_name):
    #saves the names of genomes in a data frame
    genome_names = pd.DataFrame(synbio_bacteria.index)
    # implements EC filter on the big dataset. Multiplies all of the columns by the EC weight
    print("Begining EC analysis")
    print(pref_set)
    if pref_set == 0:
        ec_filter = pd.read_csv(ec_name, delimiter = "\t", header = 0, index_col = 0)
        ec_filter.pivot(columns=["EC Number", "Scoring Based on Function"])
        print(ec_filter)
        print(ec_filter.shape)
        synbio_bacteria = synbio_bacteria.mul(ec_filter[2,:])
    ## Calculate the distances of the datafile using the pdist funciton from scipy. The intial return will only
    ## provide the upper half of the comparisons (row-wise), to crease a symetrical matrix, then create squareform
    distances_parallel = pairwise_distances(X = synbio_bacteria, metric = 'euclidean', n_jobs = 18)
    ## NaNs at this stage may indicate header name mismatch - check if EC numbers are aligning
    #
    ## If wanting to save the full distance matrix, turn the following flag on. Note: Well above 10 GBs
    #np.savetxt("Distance_matrix_Combined.txt",distances_parallel)
    print(type(distances_parallel))
    # Converts the distance matrix of synbio as a data frame. Concat binds the row names data frame with the synbio distance
    # matrix. Result should be a [2,#of total genomes]
    distances_synbio = pd.concat([genome_names, pd.DataFrame(distances_parallel[-1, :])], axis = 1)
    distances_synbio.columns = ['Genome Name', 'Calculated Distance']
    return distances_synbio

def ec_filter_preference(desired_loc):
    print("Would you like to use/upload EC filter: (Y/N) ")
    preference = input()
    if preference=="Y":
        print("Would you like to upload your own filter? ")
        personal = input()
        if personal == "Y":
            print("Enter the name of your document: " )
            filter_used = input()
            personal_filter_path = os.path.abspath(filter_used)
            if personal_filter_path != desired_loc:
                os.shutil(filter_used, desired_loc)
                default = 0
        elif personal == "N":
            print("You have selected to use the default filter: 2022_5_24_Proposed_Filter.csv")
            filter_used = "2022_5_24_Proposed_Filter.csv"
            default  = 0
        else:
            print("Error, try again!")
            filter_used = "ERROR"
    else:
        print("No EC filter has been selected. All EC numbers will have equal weight")
        default = 1
        filter_used = " "
    return default, filter_used

def pass_to_distance(big_matrix_name, genome_ID, location):
    ## Grabbing the big matrix name and the genome ID that the user provides.
    ## Returns a list of distances found
    compare_mat = read_in_binary_matrix(big_matrix_name)
    ec_pref = ec_filter_preference(location)
    setting = ec_pref[0]
    filter_name = ec_pref[1]
    distances_mat = calculate_dsitance(compare_mat, setting, filter_name)

    return distances_mat


