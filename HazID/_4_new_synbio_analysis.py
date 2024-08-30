# User input. Choose synbio file, for single synbio file
import subprocess #as subprocess; print('subprocess version:', subprocess.__version__)
import os #; print('os version:', os.__version__)
import shutil #; print('shutil version:', shutil.__version__)
import numpy as np#; print('numpy version:', np.__version__)
import re#; print('re version:', re.__version__)
import pandas as pd#; print('pandas version:', pd.__version__)
from _2_Genomic_Library_Functions import diamond_impl, genome_extractor
from _5_genomic_comparison import pass_to_distance
import sys
#from _1_competitorFind import ec_locator
##====================================================================================================================##
##This is all directory creation and file management 

sb_filename = '/projects/jodo9280/EcoDr/TestCases/GCF_002926195.1_ASM292619v1_protein.faa' #Takes in the synbio FASTA file 

file_loc = os.path.abspath(sb_filename) 

sb_name = 'New_Synbio_Analysis_Output_Binary'  


current_directory = os.getcwd()

desired_location='/projects/jodo9280/EcoDr/TestCases' 

os.chdir(desired_location) #changes current working directory to the one specific above
if os.path.exists(sb_name):
    shutil.rmtree(sb_name)
    os.mkdir(sb_name)
else:#makes a new directory called sb_name
    os.mkdir(sb_name)
desired_location = desired_location + "/" + sb_name

shutil.copy(file_loc, desired_location) #moves synbio file to Test Cases folder
matches = sb_filename + "_matches"
os.chdir(desired_location)
##====================================================================================================================##
# Calls on function from genomic data download, returns the location of the diamond results
# Calls on script _2_genomic_data_download.py
diamond_results_loc = diamond_impl(desired_location, sb_name) #runs a diamond analysis and outputs a folder location
print(diamond_results_loc)
print("DIAMOND Results Completed")
##====================================================================================================================##
# Summary Matrix Rendering-- Completes the binary matrix using the diamond hits outout (matches)
# Uses function from genomic_summary, requires the EC number full pathway which is subjected to change.
analysis_output = genome_extractor(diamond_results_loc, sb_name)
print('EC Binary scoring is complete') 
print('Analysis Output is in ', analysis_output[0])
synbio_binary = analysis_output[0]
print('EC_Binary is in ', synbio_binary)
##====================================================================================================================##
# Passes the filename of the binary matrix and the ID of organism
# Returns the distance matrix of the combined EC space for Bacteria, Archaea, and the synbio organism
# Calls on script _5_synbio_distance_matrix.py
[distance_list_for_synbio, new_loc ]= pass_to_distance(synbio_binary, sb_name, desired_location)
#synbio_binary = your binary matrix
#sb_name = 'new synbio analysis output
#desired_location = '/projects/jodo9280/EcoDr/TestCases/new symbio analysis output' 
print(distance_list_for_synbio)  # Prints distance matrix for the synbio genome
print('Synbio analysis is complete')
##====================================================================================================================##
