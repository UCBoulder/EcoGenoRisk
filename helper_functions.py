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


def input_domain():
    domain = input("Enter domain:")
    lower_case = domain.lower()
    ncbi_url=''
    if lower_case == 'archaea':
        ncbi_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt'

    elif lower_case == 'bacteria':
        ncbi_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt'

    elif lower_case == 'fungi':
        ncbi_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt'

    else:
        print("Invalid Input")

    return ncbi_url, lower_case


# Downloading/checking for assembly file
def checking_assembly_file(text, link):
    ua = UserAgent()
    header = {'User-Agent': str(ua.chrome)}
    if os.path.exists(text):
        return "File exists"
    else:
        temp_genome_list = requests.get(link, headers=header)
        temp_genome_list = requests.get(link)
        print(temp_genome_list.text)
        with open(text, 'w') as genome_list_out:
            genome_list_out.write(temp_genome_list.text)
        genome_list_out.close()


# looping through assembly summary file, finding the samples that have a complete genome, and downloading the FTP file
def file_extraction(text, dom):
    with open(text, 'r') as assembly_summary:
        assembly_summary.readline()
        for line in assembly_summary:
            output = line.split("\t")
            if re.match(output[11], "Complete Genome"):
                link = output[19]
                org_name = link
                species_name = re.search(r'(.*)/(.*)', org_name).group(2)
                url_new = link + '/' + species_name + '_protein.faa.gz'
                # wget.download(url_new)
                download_name = species_name + '_protein.faa.gz'
                destination = os.path.abspath(dom+ '_protein_file')
                source = os.path.abspath(download_name)
                # shutil.move(source, destination)
    return destination
    assembly_summary.close()


# Extracting all files and deleting the zip folder
def file_management(dest):
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
def diamond_impl(dest):
    uniprot = open('/home/anna/PycharmProjects/pythonProject/uniprot.fasta')
    reference = uniprot.read()
    if os.path.abspath('home/anna/PycharmProjects/pythonProject/DIAMOND_matches') == 'EMPTY':
        directory = '/home/anna/PycharmProjects/pythonProject/DIAMOND_matches'
        os.makedirs(directory)
    # yeet = Popen(['diamond','version'])
    # print('stdout:', yeet.stdout)
    makedb = ['diamond', 'makedb', '--', 'in', reference, '-d', 'reference']
    lib = Popen(makedb, stdin=PIPE, stdout=PIPE)
    # lib=subprocess.call(makedb)
    lib.communicate()
    counter = 0
    for item in os.listdir(dest):
        if item.endswith('.faa'):
            print(counter)
            file_path = os.path.abspath(item)
            file_name = (os.path.basename(file_path)).rsplit('.', 1)[0]
            matches = file_name + "_matches"
            # run_diamond(file_name, reference, "blastp")
            time.sleep(10)
            blastp = ['diamond', 'blastp', '-d', 'reference', '-q', file_name, '-o', matches]
            time.sleep(10)
            # subprocess.run(str(blastp), capture_output=True, text=True)
            # search=Popen(blastp, stdin=PIPE,stdout=PIPE)
            # search=subprocess.call(blastp)
            # search.communicate()
            shutil.move(os.path.abspath(matches), directory)
            counter += 1
    print(counter)
    uniprot.close()
