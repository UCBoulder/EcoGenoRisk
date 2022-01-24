from genomic_data_download import input_domain, checking_assembly_file, file_extraction, file_management, diamond_impl

path_folder = '/home/anna/PycharmProjects/pythonProject'
text_file = 'assembly_summary.txt'
domain = input_domain()
url = domain[0]
domain = domain[1]
checking_assembly_file(text_file, url)
destination= file_extraction(text_file, domain)
print(destination)

file_management(destination)
diamond_impl(destination)



print("That's all folks")
