# Creates initial files needed for analysis. This program allows you to download FASTA protein files of all complete
# organisms in a given domain, completes DIAMOND analysis on them by using the Uniprot EC number catalogue, and creates
# a binary summary matrix that shows which EC groups a certain organism has

from genomic_data_download import input_domain, checking_assembly_file, file_extraction, file_management, diamond_impl
from genomic_summary import genome_extractor, distance_based_tree

path_folder = '/home/anna/PycharmProjects/pythonProject'
text_file = 'assembly_summary.txt'
domain = input_domain()
url = domain[0]
domain = domain[1]
checking_assembly_file(text_file, url)
destination = file_extraction(text_file, domain)
file_management(destination)
new_folder = diamond_impl(destination)
output = genome_extractor(new_folder, domain)
ec_profile_matrix=output[0]
saved_file_name=output[1]
distance_based_tree(ec_profile_matrix, saved_file_name, domain)

print("That's all folks")
