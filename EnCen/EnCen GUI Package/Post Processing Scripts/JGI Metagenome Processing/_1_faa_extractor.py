###_________________________________________________________--
#
# Extracts all files from .tar compressed files. /
# parent directory is the folder the extracted file is in 
# This script is to be used before directory_clear.py


import os
import tarfile


def tar_extraction(parent_folder):
    # References those files that are compressed
    extension = ".tar.gz"
    os.chdir(parent_folder)
    # Iterates through the protein file folder, destination where the .faa files were downloaded
    for item in os.listdir(parent_folder):
        filepath = os.path.join(parent_folder, item)
        # Check for ".gz" extension, does not apply to those files that were already extracted
        if item.endswith(extension):
            print('File Found, Beginning Extraction')
            print(item)
            with tarfile.open(filepath, 'r:gz') as tar:
            # Extract all contents of the archive
                tar.extractall(path=parent_folder)
            os.remove(filepath)

def delete_files_not_ending_with(directory, extension):
    for filename in os.listdir(directory):
        if not filename.endswith(extension):
            try:
                filepath = os.path.join(directory, filename)
                os.remove(filepath)
                print('Deleted :', filepath)
            except:
                print('Skipped over', filename)


# parent_folder = '/home/anna/Documents/EnCen_JGI_Genomes'
parent_folder = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/Model Validation/Soil Rhizosphere Metagenome/img_data_295860-23'
file_extension_to_keep1 = '.faa'
tar_extraction(parent_folder)
delete_files_not_ending_with(parent_folder, file_extension_to_keep1)