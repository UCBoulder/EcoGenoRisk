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
    path = input("Enter desired pathway to make a new folder in: ")
    # os.chdir(path)
    # os.makedirs(domain)
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
    pathway = path + '/' + domain
    os.chdir(pathway)
    print(pathway)
    return [ncbi_url, lower_case]  # returns domain name and url


# Downloading/checking for assembly file
def checking_assembly_file(text, link):  # asks for assembly_summary text file name, requires download link
    print(os.getcwd())
    ua = UserAgent()
    header = {'User-Agent': str(ua.chrome)}
    if os.path.exists(text):  # checks if file exists already in order to not spam NCBI
        time.sleep(4)
        return "File exists"
    else:
        time.sleep(4)
        temp_genome_list = requests.get(link, headers=header)  # permissions for opening file through NCBI
        # temp_genome_list = requests.get(link)
        print(temp_genome_list.text)

        with open(text, 'w') as genome_list_out:
            genome_list_out.write(
                temp_genome_list.text)  # copies content from the website, pastes it into assembly_summary file
        genome_list_out.close()
        print("checking_assembly_file--success")


# looping through assembly summary file, finding the samples that have a complete genome, and downloading the FTP file
def file_extraction(text, dom):
    print(os.getcwd())
    if os.path.exists(os.getcwd() + "/" + dom + "_protein_file"):
        print('')
        destination = (os.getcwd() + "/" + dom + "_protein_file")
    else:
        os.makedirs(dom + "_protein_file")
        destination = os.path.abspath(dom + '_protein_file')  # creates a domain specified folder pathway
        with open(text, 'r') as assembly_summary:
            assembly_summary.readline()  # reads line by line, iterating through the doc
            for line in assembly_summary:
                output = line.split("\t")  # separating line by tab
                if re.match(output[11], "Complete Genome"):
                    link = output[19]  # finds ftp path
                    org_name = link
                    species_name = re.search(r'(.*)/(.*)', org_name).group(2)
                    url_new = link + '/' + species_name + '_protein.faa.gz'  # finds exact address to download the exact url of the FTP file
                    if os.path.exists(os.getcwd() + "/" + dom + "_protein_file" + "/" + species_name + "_protein.faa"):
                        print("")
                    else:
                        wget.download(url_new)
                        download_name = species_name + '_protein.faa.gz'
                        source = os.path.abspath(download_name)
                        shutil.move(source, destination)   # moves file into domain specified folder
        assembly_summary.close()
        print("file_extraction--success")
    return destination

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

    print("file_management-- success")


# DIAMOND IMPLEMENTATION
def diamond_impl(dest):  # takes domain folder path

    current_dir = os.getcwd()
    directory = current_dir.rsplit('/', 1)[0]
    print(directory)
    os.chdir(directory)
    new_destination = directory + '/DIAMOND_matches'
    print(new_destination)
    if os.path.exists(new_destination):
        print(" ")
        os.chdir(directory+'/DIAMOND_matches')
    else:
        os.makedirs('DIAMOND_matches')  # if library is present, the directory is changed to the folders path

    if os.path.exists(
            directory + '/reference.dmnd'):  # if there currently is no reference library (.dmnd), makedb creates a DIAMOND library
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
            # print(file_path)
            file_name = (os.path.basename(file_path)).rsplit('.', 1)[0]  # name of sample
            matches = file_name + "_matches"
            print(matches)
            if os.path.exists(new_destination+"/"+matches):
                print(" ")
            else:
                blastp = ['diamond', 'blastp', '-d', directory + '/reference.dmnd', '-q', file_path, '-o', matches,
                          '--max-target-seqs', '1', '--outfmt', '6']  # completes DIAMOND search by using full path
                subprocess.run(blastp)
                shutil.move(os.path.abspath(matches), new_destination)  # moves to DIAMOND folder
    print("diamond_impl--success")
    return directory
