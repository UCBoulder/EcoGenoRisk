from genomic_data_download import input_domain, checking_assembly_file, file_extraction, file_management, diamond_impl
from genomic_summary import genome_extractor

path_folder = '/home/anna/PycharmProjects/pythonProject'
text_file = 'assembly_summary.txt'
domain = input_domain()
url = domain[0]
domain = domain[1]
checking_assembly_file(text_file, url)
destination = file_extraction(text_file, domain)
file_management(destination)
ye = diamond_impl(destination)
genome_extractor(ye)
print("That's all folks")
