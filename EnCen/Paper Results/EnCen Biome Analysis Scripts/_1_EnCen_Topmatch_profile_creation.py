import gzip#; print('gzip version: ', gzip.__version__)
import os#; print('os version: ', os.__version__)
import re#; print('re version: ', re.__version__)
import shutil#; print('shutil version: ', shutil.__version__)
import requests#; print('requests version: ', requests.__version__)
import subprocess#; print('subprocess version: ', subprocess.__version__)
import time#; print('time version: ', time.__version__)
from datetime import datetime
from subprocess import PIPE, Popen
from fake_useragent import UserAgent
import os.path
from os import path
from Bio.ExPASy import Enzyme
import numpy as np 
from sklearn.metrics import pairwise_distances
import pandas as pd 
import numpy as np 



def EC_extract():
    ec_library = 'EC_library.csv' 
    ua = UserAgent()
    header = {'User-Agent': str(ua.chrome)}
    ec_url = 'https://ftp.expasy.org/databases/enzyme/enzyme.dat'
    time.sleep(4)
    ec = requests.get(ec_url, headers=header)
    if ec.status_code == 200:
        with open(ec_library, 'w+', newline='\n') as ec_file:
            ec_file.write(ec.text)
###This step creates the list 
    handle = open(ec_library)
    records = Enzyme.parse(handle)
    ecnumbers = [record["ID"] for record in records]
    #print(type(ecnumbers)) #This is a list at this point in the code
    path = os.path.abspath(ec_library)
    with open(path, 'w+', newline='\n') as csv_file:
        #csv_file = csv.writer(csv_file)
        for item in ecnumbers:
            csv_file.write(item +'\n')
    print('EC list Has Been Created')
def tsv_to_fasta():
    reference_library = 'uniprot.tsv' 
    #this creates a file named 'uniprot.tsv' on the working directory. Should show within the EcoDr/EcoDr folder. 
    ua = UserAgent()
    header = {'User-Agent': str(ua.chrome)}
    uniprot_url = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cec%2Csequence&format=tsv&query=%28%28ec%3A*%29%29+AND+%28reviewed%3Atrue%29'
    time.sleep(4)
    uniprot = requests.get(uniprot_url, headers=header)
    if uniprot.status_code == 200:
        with open(reference_library, 'w+') as reference_library:
            reference_library.write(uniprot.text)
    print('Uniprot Reference File Has Been Created')
####this step takes it the initial tsv and converts it to FASTA
    input = os.path.abspath('uniprot.tsv')
    output = 'uniprot.fasta'
    with open(input, 'r') as input_file, open(output, 'w') as output_file:
        header=next(input_file)
        for line in input_file:
            carrot = f'>{line}'
            new_row = re.sub(r'(.*?\t.*?)\t', r'\1\n', carrot, 1)
            new_row_with_tab_replaced_by_question_mark = new_row.replace("\t", "?")
            # Be careful with the symbols used for replacement, as they have to align with the genomic summary string parsing
            output_row_with_space_replaced_with_ampersand = new_row_with_tab_replaced_by_question_mark.replace("; ", ";_")
            #new_row2 = re.sub(r'(.*?\t.*?)\t', r'$', new_row, 0)
            #this is regex, for the record. 
            output_file.write(output_row_with_space_replaced_with_ampersand)
    print('FASTA has been created from TSV and is named', output)

#__________________________________________________________________________________________________________
def diamond_impl(dest, name):
    print(os.getcwd())
    matches = ''
    output_folder = dest 
    final_folder = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/Paper Results'
    print("DIAMOND library is located in: ", output_folder)
    if os.path.isfile('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/Paper Results/Uniprot_Reference_Library.dmnd'):
        print("Library Detected")
    # If not present, then creates a DIAMOND library by referencing the exact location where the Uniprot library is saved
    # If there currently is no reference library (.dmnd), then command makedb creates a DIAMOND library
    else:
        print("Creation of DIAMOND-formatted library...")
        makedb = ['diamond', 'makedb', '--in', '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/Paper Results/uniprot.fasta', '-d',
                  'Uniprot_Reference_Library.dmnd']  # Reference library full pathway
        #This is a list for the DIAMOND specific makedb function. 
        subprocess.run(makedb)
        #This simply runs the function makedb
        print("Library complete")
#^this chunk is good to go. Just makes the reference database from the uniprot fasta. 
##_______________________________________________________________# This portion does the matching. The above portion creates the reference library from the unitprot fasta
    for item in os.listdir(dest):
        # Checks for file extension
        if item.endswith('.faa') and not os.path.isfile(
                os.path.basename(os.path.abspath(item)).rsplit('.', 1)[0] + "_matches.tsv"):
            # Finds path of file
            file_path = os.path.abspath(item)
            # Finds the GCF/ASM name of the file by looking at the first part of the name before the .faa notation
            if name == "":
                print(os.path.basename(file_path).rsplit('.',1))
                file_name = (os.path.basename(file_path)).rsplit('.', 1)[0]
            else:
                file_name = name
            # New filename that ends with matches
            matches = file_name + "_matches.tsv"
            # print(matches)
            # If genome has not already undergone DIAMOND search and is currently located in the correct folder, then
            # the subprocess function will run the diamond search
            if not os.path.isfile(dest + "/" + matches) and os.path.abspath(matches) != output_folder:
                print("Processing ", file_name)
                # DIAMOND search using the full pathway of the protein files, max target sequence outputs only one best
                # match with highest e-value which represent the chance of obtaining a better random match in the same database (Buchfink et al, 2021)
                blastp = ['diamond', 'blastp', '-d', 'Uniprot_Reference_Library.dmnd', '-q', file_path, '-o', matches, 
                          '--max-target-seqs', '1', '--outfmt', '6']
                time.sleep(4)
                subprocess.run(blastp)
    

    # Moves all DIAMOND search outputs into the folder
        if item.endswith('.tsv'):
            if os.path.exists(os.path.join(final_folder, item)):
                print(f"Overwriting: {item}")
                os.remove(os.path.join(final_folder, item))
            shutil.move(os.path.abspath(item), final_folder)
    print("diamond_impl--success")
    # Returns the location of the DIAMOND matches folder
    return output_folder
#__________________________________

def genome_extractor(diamond_folder, name):
    os.chdir(diamond_folder)
    print(os.getcwd())
    # Opens the list of of EC numbers
    ec_open = np.loadtxt('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/Paper Results/EC_library.csv',
                         dtype='str')
    big_matrix = ["Name_of_Genome"]
    file_name = name + '_functional_profile'
    new_dir = diamond_folder + '/' + file_name
    # Checks to see if the document already exists using full pathway name
    if os.path.exists(new_dir):
        pass
        #print("Summary Matrix exists")
        #return [new_dir, file_name]
    else:
        for ec_force in ec_open:
            # Creates a horizontal header of all of the EC names
            big_matrix.append(ec_force)
   
        for item in os.listdir(diamond_folder):
            if item.endswith("_matches.tsv"):
                print(item)
                genome = [item] #Turns the GCF's into a list, where the GCF names in the matrix come from. 
                genome_runner_ec = [item] #Turns the GCF's into a list, where the EC is appended
                # Iterates through all of the EC numbers (1:8197)
                
                GCF = open(item, 'r') # CBM Added
                
                for line in GCF: # CBM Added
                    no_tab = line.split('\t')
                    first_ec = no_tab[1].split("?")
                    separate_ec = first_ec[1].split(";_")
                    genome_runner_ec.append(separate_ec)
                    print("Appending...")

                for ec in ec_open:
                 
                    ec_now = 0
                    if [ec] in genome_runner_ec:
                        ec_now = 1

                    # 1 or 0 will be appended to the summary matrix for each EC value in the list
                    genome.append(ec_now)
                    #print(genome)
                # Vertical stacking occurs for each genome in the DIAMOND output folder
                big_matrix = np.vstack([big_matrix, genome])
        #print(big_matrix)
        # Saves matrix as a text file for further analysis
        np.savetxt(file_name, big_matrix, fmt='%s')
        # Returns the location of the summary matrix and the name of the file
        # if not os.path.exists(os.path.abspath(file_name)):
        #     shutil.move(os.path.abspath(file_name),'/home/anna/Documents/UBA6164')
        # else:
        #     print('')
        print(new_dir)
    return [new_dir, file_name]
# ______________________________________________________________________________________
synbio = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/Paper Results/Biome Analysis Results/Lake_Sediment'
name = '3300045561_8(Lake_sediment_Top_Match)'
syn_folder_name = '3300045561_8(Lake_sediment_Top_Match)'
# desired_location2 = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/Paper Results/Biome Analysis Results/Activated Sludge/New_results'

os.chdir(synbio)

# EC_extract()
# tsv_to_fasta()
 
diamond_syn = diamond_impl(synbio, name)

output2 = genome_extractor(diamond_syn, name)
