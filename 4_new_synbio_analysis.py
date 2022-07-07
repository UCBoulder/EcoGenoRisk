# User input. Choose synbio file, for single synbio file
import subprocess, import os, import shutil, import numpy as np, import re, import pandas as pd
from 2_genomic_data_download_editted import diamond_impl
from 3_genomic_summary import genome_extractor
from 5_synbio_distance_matrix import pass_to_distance
#Turn off when testing
print("Welcome to EcoGeno! Enter the name of your FASTA file: ")
sb_filename = input()
print("Enter where you want results saved: ")
desired_location = input()
print("Enter organism ID, no punctuation: ")
sb_name = input()
desired_location = "/home/anna/PycharmProjects/pythonProject/Bacteria"
#Turn on when testing
# sb_name = "Hello"
# sb_filename = "Test.faa"
matches = sb_filename+"_matches"
#Calls on function from genomic data download, returns the location of the diamond results
diamond_results_loc = diamond_impl(desired_location)
print("DIAMOND Results Completed")
## Summary Matrix Rendering-- Completes the binary matrix using the diamond hits outout (matches)
# Uses function from genomic_summary, requires the EC number full pathway which is subjected to change.
analysis_output = genome_extractor(diamond_results_loc)
name = analysis_output[1]
print("EC Binary Scoring is Complete")
# Passes the filename of the binary matrix and the ID of organism
# Returns the distance matrix of the combined ec space for bacteria, archaea, and synbio organism
distance_list_for_synbio = pass_to_distance(name, sb_name, desired_location)
print(distance_list_for_synbio)             #Prints distance matrix for the synbio genome
np.savetxt("DistanceListforSynbio_Class.txt", distance_list_for_synbio, fmt='%s', delimiter= '\t')   #saves file of the synbio distance matrix
##Further processing can be found in