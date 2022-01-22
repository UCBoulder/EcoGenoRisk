from helper_functions import input_domain, checking_assembly_file, file_extraction, file_management, diamond_impl

path_folder = '/home/anna/PycharmProjects/pythonProject'
text_file = 'assembly_summary.txt'
input_domain()
url = input_domain()[0]
domain = input_domain()[1]
checking_assembly_file(text_file, url)
destination= file_extraction(text_file, domain)
file_management(destination)
# diamond_impl(destination)

print("That's all folks")
