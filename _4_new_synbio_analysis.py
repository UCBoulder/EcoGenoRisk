# User input. Choose synbio file, for single synbio file
import subprocess #as subprocess; print('subprocess version:', subprocess.__version__)
import os #; print('os version:', os.__version__)
import shutil #; print('shutil version:', shutil.__version__)
import numpy as np#; print('numpy version:', np.__version__)
import re#; print('re version:', re.__version__)
import pandas as pd#; print('pandas version:', pd.__version__)
from _2_genomic_data_download import diamond_impl
from _3_genomic_summary import genome_extractor
from _5_synbio_distance_matrix import pass_to_distance
from _1_competitorFind import ec_locator
##====================================================================================================================##
# Asks for user preferences and synbio-protein file related info
print("Welcome to EcoGeno! Enter the name of your FASTA file: ")
sb_filename = input()
file_loc = os.path.abspath(sb_filename)
print("Enter organism ID, no punctuation (no space in between words): ")
sb_name = input()
print("Enter where you want results saved: ")
os.chdir(desired_location)
os.makedirs(sb_name)
desired_location = desired_location + "/" + sb_name
print(desired_location)
shutil.move(file_loc, desired_location)
matches = sb_filename + "_matches"
os.chdir(desired_location)
##====================================================================================================================##
# Calls on function from genomic data download, returns the location of the diamond results
# Calls on script _2_genomic_data_download.py
diamond_results_loc = diamond_impl(desired_location, sb_name)
print(diamond_results_loc)
print("DIAMOND Results Completed")
##====================================================================================================================##
# Summary Matrix Rendering-- Completes the binary matrix using the diamond hits outout (matches)
# Uses function from genomic_summary, requires the EC number full pathway which is subjected to change.
# Calls on script _3_genomic_summary.py
analysis_output = genome_extractor(diamond_results_loc, sb_name)
print(analysis_output)
ec_binary = analysis_output[0]
print("EC Binary Scoring is Complete")
##====================================================================================================================##
# Passes the filename of the binary matrix and the ID of organism
# Returns the distance matrix of the combined EC space for Bacteria, Archaea, and the synbio organism
# Calls on script _5_synbio_distance_matrix.py
[distance_list_for_synbio, new_loc ]= pass_to_distance(ec_binary, sb_name, desired_location)
print(distance_list_for_synbio)  # Prints distance matrix for the synbio genome
print('Synbio analysis is complete')
##====================================================================================================================##
# Sends synbio name and pathway of saved files to CompetitorFind scripts
print = new_loc
ec_locator(sb_name, new_loc)