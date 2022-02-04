import gzip
import os
import re
import shutil
import requests
import wget
import subprocess
from subprocess import PIPE, Popen
import time
from fake_useragent import UserAgent


# Asks for input for domain, returns a specific url
def input_domain():
    domain = input("Enter domain:")
    lower_case = domain.lower()
    ncbi_url = ''
    if lower_case == 'archaea':
        ncbi_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt'

    elif lower_case == 'bacteria':
        ncbi_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt'

    elif lower_case == 'fungi':
        ncbi_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt'

    else:
        print("Invalid Input")

    return ncbi_url, lower_case  # returns domain name and url


# Downloading/checking for assembly file
def checking_assembly_file(text, link):  # asks for assembly_summary text file name, requires download link
    ua = UserAgent()
    header = {'User-Agent': str(ua.chrome)}
    if os.path.exists(text):  # checks if file exists already in order to not spam NCBI
        return "File exists"
    else:
        temp_genome_list = requests.get(link, headers=header)  # permissions for opening file through NCBI
        # temp_genome_list = requests.get(link)
        print(temp_genome_list.text)
        with open(text, 'w') as genome_list_out:
            genome_list_out.write(
                temp_genome_list.text)  # copies content from the website, pastes it into assembly_summary file
        genome_list_out.close()


# looping through assembly summary file, finding the samples that have a complete genome, and downloading the FTP file
def file_extraction(text, dom):
    with open(text, 'r') as assembly_summary:
        assembly_summary.readline()  # reads line by line, iterating through the doc
        for line in assembly_summary:
            output = line.split("\t")  # separating line by tab
            if re.match(output[11], "Complete Genome"):
                link = output[19]  # finds ftp path
                org_name = link
                species_name = re.search(r'(.*)/(.*)', org_name).group(2)
                url_new = link + '/' + species_name + '_protein.faa.gz'  # finds exact address to download the exact url of the FTP file
                # wget.download(url_new)
                download_name = species_name + '_protein.faa.gz'
                destination = os.path.abspath(dom + '_protein_file')  # creates a domain specified folder pathway
                source = os.path.abspath(download_name)
                # shutil.move(source, destination)   # moves file into domain specified folder
    return destination
    assembly_summary.close()


# Extracting all files and deleting the zip folder
def file_management(dest):  # takes the destination of the domain folder
    dir_name = 'x'
    extension = ".gz"
    file_name = ''
    os.chdir(dest)
    for item in os.listdir(dest):  # loop through items in dir
        if item.endswith(extension):  # check for ".gz" extension
            gz_name = os.path.abspath(item)  # get full path of files
            file_name = (os.path.basename(gz_name)).rsplit('.', 1)[0]  # get file name for file within
            with gzip.open(gz_name, "rb") as f_in, open(file_name, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(gz_name)  # delete zipped file


# DIAMOND IMPLEMENTATION
def diamond_impl(dest):  # takes domain folder path

    directory='/home/anna/PycharmProjects/pythonProject/DIAMOND_matches'
    if os.path.exists(directory):
        print(directory)
        os.chdir(directory)
    else:
        os.chdir('/home/anna/PycharmProjects/pythonProject')
        os.makedirs('DIAMOND_matches')  # if library is present, the directory is changed to the folders path
        os.chdir(directory)
        print(directory)

    if os.path.exists('/home/anna/PycharmProjects/pythonProject/DIAMOND_matches/reference.dmnd'):  # if there currently is no reference library (.dmnd), makedb creates a DIAMOND library
        print("Library Detected")
    else:
        makedb = ['diamond', 'makedb', '--in', '/home/anna/PycharmProjects/pythonProject/uniprot.fasta', '-d',
                  'reference']  # reference library full pathway
        subprocess.run(makedb)
        print("Library complete")

    os.chdir(dest)
    for item in os.listdir(dest):  # for item in domain_folder, if item is a faa file, complete diamond analysis
        if item.endswith('.faa'):
            file_path = os.path.abspath(item)
            print(file_path)
            file_name = (os.path.basename(file_path)).rsplit('.', 1)[0]  # name of sample
            matches = file_name + "_matches"
            if os.path.exists(directory+"/"+matches):
                print(" ")
            else:
                blastp = ['diamond', 'blastp', '-d', '/home/anna/PycharmProjects/pythonProject/DIAMOND_matches/reference.dmnd', '-q', file_path, '-o', matches, '--max-target-seqs', '1', '--outfmt', '6']  # completes DIAMOND search by using full path
                subprocess.run(blastp)
                shutil.move(os.path.abspath(matches), directory)  # moves to DIAMOND folder
    return directory
