# User input. Choose synbio file, for single synbio file
import subprocess; import os; import shutil; import numpy as np; import re; import pandas as pd
from _2_genomic_data_download_editted import diamond_impl
from _3_genomic_summary import genome_extractor
from _5_synbio_distance_matrix import pass_to_distance
#Turn off when testing
# print("Welcome to EcoGeno! Enter the name of your FASTA file: ")
# sb_filename = input()
# print("Enter organism ID, no punctuation: ")
# sb_name = input()
# print("Enter where you want results saved: ")
# desired_location = input()
# os.chdir(desired_location)
# os.makedirs(sb_name)
# desired_location = desired_location+"/"+sb_name
# shutil.move(os.path.abspath(sb_filename), desired_location)
#Turn on when testing
desired_location = "/home/anna/PycharmProjects/pythonProject/20220710_testing"
#os.chdir(desired_location)
sb_name = "Synbio1"
#sb_filename = "GCF_000009545.1_ASM954v1_protein.faa"
sb_filename = 'GCF_000012805.1_ASM1280v1_protein.faa'
# Testing of different files for distance matrix output
#sb_name = "Desulfurococcus"
#sb_filename = "GCF_000186365.1_ASM18636v1_protein.faa"
#sb_name = "Neisseria"
#sb_filename = "GCF_000011705.1_ASM1170v1_protein.faa"
#sb_name = "Serratia"
#sb_filename = "GCF_000214195.1_ASM21419v1_protein.faa"
# Issues with file management and moving to appropriate places, construct code after the distance matrix comes out right
#desired_location = "/home/anna/PycharmProjects/pythonProject/"+sb_filename
# if not os.path.isdir(desired_location):
#     os.makedirs(sb_filename)
#     shutil.move(os.path.abspath(sb_filename), desired_location)

#matches = sb_filename+"_matches"
#Calls on function from genomic data download, returns the location of the diamond results
## Calls on script 2
#diamond_results_loc = diamond_impl(desired_location, sb_name)
#print("DIAMOND Results Completed")
## Summary Matrix Rendering-- Completes the binary matrix using the diamond hits outout (matches)
# Uses function from genomic_summary, requires the EC number full pathway which is subjected to change.
## Calls on script 3
analysis_output = genome_extractor(desired_location)
print(analysis_output)
name = analysis_output[1]
print(name)
print("EC Binary Scoring is Complete")
# Passes the filename of the binary matrix and the ID of organism
# Returns the distance matrix of the combined ec space for bacteria, archaea, and synbio organism
## Calls on script 5
distance_list_for_synbio = pass_to_distance(name, sb_name, desired_location)
print(distance_list_for_synbio)             #Prints distance matrix for the synbio genome
##Further processing can be found in