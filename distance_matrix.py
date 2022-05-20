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
    synbio_bacteria = bacteria_binary.append(binary_matrix_table)       #Verically adding synbio matrix to the overall matrix in order to complete comparison
    print("Shape of combined matrix using append is: ", np.shape(synbio_bacteria))
    print("Merging of bacteria data and synbio data is complete")
    return synbio_bacteria

def calculate_dsitance(synbio_bacteria):
    ## Calculate the distances of the datafile using the pdist funciton from scipy. The intial return will only
    ## provide the upper half of the comparisons (row-wise), to crease a symetrical matrix, then create squareform
    distances_parallel = pairwise_distances(X = synbio_bacteria, metric = 'euclidean', n_jobs = 18)     #Outputs the distances found, Cannot complete analysis since input contains NaN
    np.savetxt("Distance_matrix_Combined.txt",distances_parallel)
    return distances_parallel

#np.savetxt("Genome_to_genome_difference_comparison.txt", distances_parallel, fmt='%s')

#read_in_binary_matrix("bacteria_big_matrix_tab.txt")

def pass_to_distance(big_matrix_name, genome_ID):
    ## Grabbing the big matrix name and the genome ID that the user provides.
    ## Returns a list of distances found

    compare_mat = read_in_binary_matrix(big_matrix_name)
    distances_mat = calculate_dsitance(compare_mat)
    return distances_mat[-1,:]
