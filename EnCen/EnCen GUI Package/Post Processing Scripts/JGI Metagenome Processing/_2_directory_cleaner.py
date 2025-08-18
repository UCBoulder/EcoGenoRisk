###This code iterates through a directory and deletes all files that are not a .faa
import os

def delete_files_not_ending_with(directory, extension):
    for filename in os.listdir(directory):
        if not filename.endswith(extension):
            try:
                filepath = os.path.join(directory, filename)
                os.remove(filepath)
                print('Deleted :', filepath)
            except:
                print('Skipped over', filename)
# Provide the directory path
directory_path = '/home/anna/Documents/Soil_extracted_genomes'

# Provide the file extension you want to keep
file_extension_to_keep = '.faa'

delete_files_not_ending_with(directory_path, file_extension_to_keep)