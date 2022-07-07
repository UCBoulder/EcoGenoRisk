# Creates initial files needed for analysis. This program allows you to download FASTA protein files of all complete
# organisms in a given domain, completes DIAMOND analysis by using the Uniprot EC number catalogue, and creates
# a binary summary matrix that shows which EC groups a certain organism has
#Inputs: user-inputted type of domain
#Outputs: FASTA protein files for complete genome assembly, DIAMOND searches, EC binary summary matrix
from genomic_data_download import input_domain, checking_assembly_file, file_extraction, file_management, diamond_impl
from genomic_summary import genome_extractor, distance_based_tree

#Creates NCBI ftp domain name for assembly summary document, saves domain name for namingg convention
domain = input_domain()
url = domain[0]
new_domain = domain[1]

#Assembly summary text file will be used as a reference sheet for protein file download. Only complete genomes will be
# downloaded and used for analysis
text_file = "assembly_summary.txt"
#Checks to see if the assembly summary file exists on path, if not present, creates a new file
checking_assembly_file(text_file, url)

#Checks to see if individual organism protein files exist, if not present, downloads new file
destination = file_extraction(text_file, new_domain)

#Unzips and extracts all protein files into desired location
file_management(destination)

#Checks for DIAMOND-based library and DIAMOND matches output, if not present, runs DIAMOND
new_folder = diamond_impl(destination)

#Iterates from DIAMOND matches and extracts the EC numbers present, constructs a binary EC summary matrix
output = genome_extractor(new_folder)
ec_profile_matrix=output[0]
saved_file_name=output[1]

print("Library Implementation is Complete")
