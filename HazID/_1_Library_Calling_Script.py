# Creates initial files needed for analysis. This program allows you to download FASTA protein files of all complete
# organisms in a given domain, completes DIAMOND analysis by using the Uniprot EC number catalogue, and creates
# a binary summary matrix that shows which EC groups a certain organism has
# Inputs: user-inputted type of domain. 
# Outputs: FASTA protein files for complete genome assembly, DIAMOND searches, EC binary summary matrix
import os
from _2_Genomic_Library_Functions import input_domain, checking_assembly_file, file_extraction, file_management, tsv_to_fasta, EC_extract, diamond_impl, genome_extractor
##===================================================================================================================##
# Creates NCBI ftp domain name for assembly summary document, saves domain name for naming convention
print('Welcome to EcoGenoRisk! Before we do anything else, we need to create a Uniprot Fasta file for reference and an updated list of EC Numbers. Please Hold.')
tsv_to_fasta()
EC_extract()
out_domain = input_domain()
url = out_domain[0]
folder1 = out_domain[1]
##===================================================================================================================##
# Assembly summary text file will be used as a reference sheet for protein file download. Only complete genomes will be
# downloaded and used for analysis
text_file = "Genome_Assembly_Summary.txt"
# Checks to see if the assembly summary file exis_1ts on path, if not present, then creates a new file
checking_assembly_file(text_file, url, folder1)
##===================================================================================================================##
# Checks to see if individual organism protein files exist, if not present, then downloads new file
destination = file_extraction(text_file, folder1) #destination being FASTA_&_DIAMOND folder
##===================================================================================================================##
# Unzips and extracts all protein files into desired location
file_management(destination)
##===================================================================================================================##
# Checks for DIAMOND - based library and DIAMOND matches output.If not present, then runs DIAMOND
os.chdir(destination) #destination being FASTA_&_DIAMOND folder
diamond_folder = diamond_impl(destination, '')
##===================================================================================================================##
#Iterates from DIAMOND matches and extracts the EC numbers present, constructs a binary EC summary matrix
output = genome_extractor(destination, '')
ec_profile_matrix=output[0]
saved_file_name=output[1]
##===================================================================================================================##
print("Library Implementation is Complete")
