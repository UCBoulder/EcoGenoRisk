# User input. Choose synbio file, for single synbio file
import subprocess; import os; import shutil; import numpy as np; import re; import pandas as pd
from _2_genomic_data_download_editted import diamond_impl
from _3_genomic_summary import genome_extractor
from _5_synbio_distance_matrix import pass_to_distance
#Turn off when testing
print("Welcome to EcoGeno! Enter the name of your FASTA file: ")
sb_filename = input()
print("Enter organism ID, no punctuation: ")
sb_name = input()
print("Enter where you want results saved: ")
#desired_location = input()
#os.chdir(desired_location)
# os.makedirs(sb_name)
# desired_location = desired_location+"/"+sb_name
# shutil.move(os.path.abspath(sb_filename), desired_location)
#Turn on when testing
# Issues with file management and moving to appropriate places, construct code after the distance matrix comes out right
desired_location = "/home/anna/PycharmProjects/HazID"
# if not os.path.exists(desired_location):
#      os.makedirs(sb_name)
#      print(desired_location)
#      shutil.move(os.path.abspath(sb_filename), desired_location)
# else:
#      print('Folder Detected')
# os.chdir(desired_location)

#matches = sb_filename+"_matches"
#Calls on function from genomic data download, returns the location of the diamond results
## Calls on script 2
#diamond_results_loc = diamond_impl(desired_location, sb_name)
#print(diamond_results_loc)
#print("DIAMOND Results Completed")
## Summary Matrix Rendering-- Completes the binary matrix using the diamond hits outout (matches)
# Uses function from genomic_summary, requires the EC number full pathway which is subjected to change.
## Calls on script 3
#analysis_output = genome_extractor(diamond_results_loc, sb_name)
#print(analysis_output)
#ec_binary = analysis_output[0]
print("EC Binary Scoring is Complete")
ec_binary = '/home/anna/PycharmProjects/HazID/DIAMOND_matches/Botulism'
# Passes the filename of the binary matrix and the ID of organism
# Returns the distance matrix of the combined ec space for bacteria, archaea, and synbio organism
## Calls on script 5
distance_list_for_synbio = pass_to_distance(ec_binary, sb_name, desired_location)
print(distance_list_for_synbio)             #Prints distance matrix for the synbio genome


