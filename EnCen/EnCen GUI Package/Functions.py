# %%
# -*- coding: utf-8 -*-
#Initial Imports
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
import tkinter as tk
from tkinter import filedialog
import streamlit as st 
import ast

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

@st.cache_resource()
def diamond_impl(dest, name, reference):
    # st.write(os.getcwd())
    os.chdir(dest)
    matches = ''
    output_folder = dest 
    # final_folder = '/home/anna/Documents/JGI_soil_genomes/diamond_analysis_output'
    print("DIAMOND library is located in: ", output_folder)
    if os.path.isfile(dest + '/Uniprot_Reference_Library.dmnd'):
        st.success(":green[Diamond Reference Library Detected]")
    # If not present, then creates a DIAMOND library by referencing the exact location where the Uniprot library is saved
    # If there currently is no reference library (.dmnd), then command makedb creates a DIAMOND library
    else:
        st.write("Creation of DIAMOND-formatted library...")
        makedb = ['diamond', 'makedb', '--in', reference, '-d',
                  'Uniprot_Reference_Library.dmnd']  # Reference library full pathway
        #This is a list for the DIAMOND specific makedb function. 
        subprocess.run(makedb)
        # shutil.move(os.path.abspath('Uniprot_Reference_Library.dmnd'), dest)
        #This simply runs the function makedb
        st.success("Diamond Reference Library Generated")
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
                blastp = ['diamond', 'blastp', '--quiet', '-d',  dest + '/' + 'Uniprot_Reference_Library.dmnd', '-q', file_path, '-o', matches, 
                          '--max-target-seqs', '1', '--outfmt', '6']
                time.sleep(4)
                subprocess.run(blastp)
    

    # Moves all DIAMOND search outputs into the folder
        # if item.endswith('.tsv'):
        #     if os.path.exists(os.path.join(dest, item)):
        #         print(f"Overwriting: {item}")
        #         os.remove(os.path.join(final_folder, item))
        #     shutil.move(os.path.abspath(item), final_folder)
    print("diamond_impl--success")
    # Returns the location of the DIAMOND matches folder
    return output_folder

@st.cache_data()
def genome_extractor_ref(name, home_dir):

    os.chdir(home_dir)
    print(os.getcwd())
    # Opens the list of of EC numbers
    ec_open = np.loadtxt(home_dir + '/EC_library.csv',
                         dtype='str')
    big_matrix = ["Name_of_Genome"]
    file_name = name + '_functional_profile'
    new_dir = home_dir + '/' + file_name
    # Checks to see if the document already exists using full pathway name
    if os.path.exists(new_dir):
        pass
        #print("Summary Matrix exists")
        #return [new_dir, file_name]
    else:
        for ec_force in ec_open:
            # Creates a horizontal header of all of the EC names
            big_matrix.append(ec_force)
   
        for item in os.listdir(home_dir):
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
                    genome.append(ec_now)
                  
                big_matrix = np.vstack([big_matrix, genome])


        np.savetxt(file_name, big_matrix, fmt='%s')
      
    return [new_dir, file_name]

@st.cache_data()
def genome_extractor_syn(diamond_folder, name, home_dir):
    os.chdir(diamond_folder)
    print(os.getcwd())
    # Opens the list of of EC numbers
    ec_open = np.loadtxt(home_dir + '/EC_library.csv',
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
                
                for line in GCF: # CBM Added, this is chunk that extracts the EC's from _matches 
                    no_tab = line.split('\t')
                    first_ec = no_tab[1].split("?")
                    separate_ec = first_ec[1].split(";_")
                    genome_runner_ec.append(separate_ec)
                    print("Appending...")

                for ec in ec_open:
                 
                    ec_now = 0
                    if [ec] in genome_runner_ec:
                        ec_now = 1
                    genome.append(ec_now)
                  
                big_matrix = np.vstack([big_matrix, genome])

        np.savetxt(file_name, big_matrix, fmt='%s')
        # if not os.path.exists(os.path.abspath(file_name)):
        #     shutil.move(os.path.abspath(file_name), home_dir)
        # else:
        #     print('File already Exists')
    return [new_dir, file_name]
##_____________________________________________________________________________Dividing between functional profiles and Scoring Functions
def genome_to_genome_diffcomp(synbio_binary, domain_binary):
    names_of_orgs = pd.DataFrame(domain_binary.index) 
    diff = pd.DataFrame(abs(domain_binary.values - synbio_binary.values))
    row_sum = diff.sum(axis=1)
    df1 = pd.DataFrame(row_sum)
    difference_based_comparison = pd.concat([names_of_orgs, df1], axis=1)
    difference_based_comparison.columns = ['Metagenome Bin ID', 'Difference Score']
    difference_based_comparison = difference_based_comparison.sort_values(by='Difference Score',
                                                                          ignore_index=True).reset_index(drop=True)
    difference_based_comparison.to_csv('Absolute_Difference_Comparison_Score.txt', header=True, index=True, sep='\t')

def read_in_binary_matrix(synbio_binary, sb_name, mg_to_analyze):
    # Converts synbio summary matrix into a dataframe
    synbio_binary = pd.read_csv(synbio_binary, delim_whitespace=True, header=0)
    # print(synbio_binary)
    synbio_binary = synbio_binary.set_index('Name_of_Genome')
    #index is a label for all rows - allows the two seperate dataframes to come together since they share an index
    # print(synbio_binary)
    # print(sb_name, " size of ", np.shape(synbio_binary), " successfully imported.")
    # Opens the matrix that includes the Bacteria and Archaea summary result
    domain_binary = pd.read_csv('/home/anna/Documents/JGI_soil_genomes/functional_profiles/' + mg_to_analyze + '_metagenome_functional_profile',
                                  delim_whitespace=True, header=0)
    #this is from chunk 1, and is the EC_Binary we generated earlier
    domain_binary = domain_binary.set_index('Name_of_Genome')
    # print(domain_binary)
    # Sends to a function for direct genome to genome comparison based on EC summary matrix
    genome_to_genome_diffcomp(synbio_binary, domain_binary)
    synbio_bacteria = pd.concat([domain_binary, synbio_binary])
    # Verically adding synbio matrix to the overall matrix to complete the comparison, therefore last index
    # should represent the synbio genome that is tested and returned as distances_mat[-1,:] in pass_to_distance.
    # Note, headers are lost and need to be directly passed
    # print("Shape of combined matrix using append is: ", np.shape(synbio_bacteria))
    print("Merging of bacteria data and synbio data is complete")
    doc_name_1 = 'combined_synbio_metagenome_binary_matrix.txt'
    # Removes any duplicates
    all_matrix = synbio_bacteria[~synbio_bacteria.index.duplicated(keep='first')]
    # print(all_matrix)
    # all_matrix.to_csv(doc_name_1, header=True, index=True, sep='\t')
    return all_matrix


def diff_score(synbio_binary, domain_binary):
    # Converts synbio summary matrix into a dataframe
    synbio_binary = pd.read_csv(synbio_binary, delimiter=" ", header=0)
    # print(synbio_binary)
    synbio_binary = synbio_binary.set_index('Name_of_Genome')
    #index is a label for all rows - allows the two seperate dataframes to come together since they share an index
    # print(synbio_binary)
    # print(sb_name, " size of ", np.shape(synbio_binary), " successfully imported.")
    # Opens the matrix that includes the Bacteria and Archaea summary result
    domain_binary = pd.read_csv(domain_binary,
                                  delimiter=" ", header=0)
    #this is from chunk 1, and is the EC_Binary we generated earlier
    domain_binary = domain_binary.set_index('Name_of_Genome')
    # print(domain_binary)
    # Sends to a function for direct genome to genome comparison based on EC summary matrix
    genome_to_genome_diffcomp(synbio_binary, domain_binary)
    synbio_bacteria = pd.concat([domain_binary, synbio_binary])
    # Verically adding synbio matrix to the overall matrix to complete the comparison, therefore last index
    # should represent the synbio genome that is tested and returned as distances_mat[-1,:] in pass_to_distance.
    # Note, headers are lost and need to be directly passed
    # print("Shape of combined matrix using append is: ", np.shape(synbio_bacteria))
    print("Merging of bacteria data and synbio data is complete")
    doc_name_1 = 'combined_synbio_metagenome_binary_matrix.txt'
    # Removes any duplicates
    all_matrix = synbio_bacteria[~synbio_bacteria.index.duplicated(keep='first')]
    # print(all_matrix)
    # all_matrix.to_csv(doc_name_1, header=True, index=True, sep='\t')
    return all_matrix

def genome_to_genome_diffcomp(synbio_binary, domain_binary):
    names_of_orgs = pd.DataFrame(domain_binary.index) 
    diff = pd.DataFrame(abs(domain_binary.values - synbio_binary.values)) 
    row_sum = diff.sum(axis=1)
    df1 = pd.DataFrame(row_sum)
    difference_based_comparison = pd.concat([names_of_orgs, df1], axis=1)
    difference_based_comparison.columns = ['Organisms Compared to Synbio', 'Difference Score']
    difference_based_comparison = difference_based_comparison.sort_values(by='Difference Score',
                                                                          ignore_index=True).reset_index(drop=True)
    difference_based_comparison.to_csv('Absolute_Difference_Comparison_Score.txt', header=True, index=True, sep='\t')

 

def upload_file(home_dir, mg_to_analyze):
            root = tk.Tk()
            root.withdraw() # Hide the main window
            file_paths = filedialog.askopenfilenames(initialdir= home_dir, title = 'Select the ' + mg_to_analyze + ' metagenome .faa files you would like to analyze against')
            return file_paths

def upload_file2(home_dir, mg_to_analyze):
        root = tk.Tk()
        root.withdraw() # Hide the main window
        file_paths = filedialog.askopenfilenames(initialdir= home_dir, title = 'Select the synbio.faa files you would like to analyze against the ' + mg_to_analyze + ' metagenome')
        return file_paths

def move_files_to_folder(file_paths, destination_folder):
        for file_path in file_paths:
            if os.path.isfile(file_path):
                shutil.move(file_path, destination_folder)

def EC_corrector(home_dir):
    with open('EC_List_For_parsing.txt', 'r') as q:
        content = q.read()
        meta_EC = ast.literal_eval(content)
    for file in os.listdir(home_dir): 
        if 'functional_profile' in file:
            path = os.path.join(home_dir, file)
            synbio_df = pd.read_csv(path, delimiter=" ", header = 0)
            filtered_df = synbio_df.filter(items = meta_EC, axis = 1)
            filtered_df.to_csv('BioPesticide Parsed.txt', sep=' ', index = False)
            
  