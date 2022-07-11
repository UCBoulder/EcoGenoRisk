import gzip; import os; import re; import shutil; import requests; import wget; import subprocess; import time
from datetime import datetime
from subprocess import PIPE, Popen
from fake_useragent import UserAgent
import os.path
from os import path

# Asks for input for domain, returns a specific NCBI RefSeq URL for protein file download
# Stores the name of domain for future naming convention
def input_domain():
    #Creates a timestamp in form of date,month,year for naming convention
    str_date = date_time.strftime("%d-%m-%Y")
    domain = input("Enter domain:")
    folder1_name = str_date+"_"+domain
    path = input("Enter desired pathway to make a new folder in: ")
    #Changes path to this new location to store all of the files in
    os.chdir(path)
    # Creates a new folder in this path. Folder name will start with the timestamp and end with the domain name
    #All new files/folders will be created in this folder
    os.makedirs(folder1_name)
    #converts all domain input for consistency
    lower_case = domain.lower()
    ncbi_url = ''
    # Matches domain name, returns NCBI RefSeq ftp domain name for protein file download
    if lower_case == 'archaea':
        ncbi_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt'
    elif lower_case == 'bacteria':
        ncbi_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt'
    elif lower_case == 'fungi':
        ncbi_url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/fungi/assembly_summary.txt'
    else:
        print("Invalid Input")
    #Changes pathway to this address
    pathway = path + '/' + folder1_name
    os.chdir(pathway)
    print("All new files will be saved in: ",pathway)
    #returns the new folder name (time stamp+ domain name) and the ncbi url
    return [ncbi_url, folder1_name]

#Writes a new document and copies information from NCBI RefSeq
# Input requires assembly summary download link and textfile name for assembly summary
# Returns an assembly summary sheet that will be used as reference document for finding complete bacterial genomes
def checking_assembly_file(text, link):  # asks for assembly_summary text file name, requires download link
    # prints current directory
    print(os.getcwd())
    # Header allows for retriving and and parsing for a web browser
    # Exmaples of user agent function for web scraping: https://stackoverflow.com/questions/27652543/how-to-use-python-requests-to-fake-a-browser-visit-a-k-a-and-generate-user-agent
    ua = UserAgent()
    header = {'User-Agent': str(ua.chrome)}
    # checks if file exists already in order to not spam NCBI
    if os.path.exists(text):
        return "Assembly Summary File Already Exists"
    else:
        # pauses run for 4 seconds to avoid spam
        time.sleep(4)
        # permissions for opening file through NCBI
        temp_genome_list = requests.get(link, headers=header)
        print(temp_genome_list.text)
        #opens the empty assembly summary text file
        with open(text, 'w') as genome_list_out:
            # copies content from the website, pastes it into assembly_summary file
            genome_list_out.write(
                temp_genome_list.text)
        #closes out document
        genome_list_out.close()
        return "checking_assembly_file--success"

#Creates a protein file folder for the most recent inquiry
# Looping through assembly summary file, finding the samples that have a complete genome, and downloading the FTP file
# Inputs requires the filled out assembly summary text file (derived from checking_assembly_file) and domain name (derived from input_domain)
# Outputs the location of the FASTA protein files for all complete genome organisms
def file_extraction(text, dom):
    # prints current directory
    print(os.getcwd())
    #Checks to see if a protein file specific folder already exist for the most recent domain inquiry
    # Prints a space if the folder is already present
    if os.path.exists(os.getcwd() + "/" + dom + "_protein_file"):
        print('')
        destination = (os.getcwd() + "/" + dom + "_protein_file")
    else:
        # if folder is not present, creates a protein file specific folder and saves the pathway of folder as the
        # desired destination for protein files
        os.makedirs(dom + "_protein_file")
        destination = os.path.abspath(dom + '_protein_file')  # creates a domain specified folder pathway
        with open(text, 'r') as assembly_summary:
            # reads one line at a time, iterating through the doc
            assembly_summary.readline()
            for line in assembly_summary:
                #Splits line by tab in order to reference a specific column in the document
                output = line.split("\t")
                #References specific column for assembly level, if status is complete, extracts the ftp path
                if re.match(output[11], "Complete Genome"):
                    #Extracting the ftp path link, this website takes you to index of all genomes that are available for download
                    link = output[19]
                    org_name = link
                    # Species name is found by extracting the genome accession number/ASM name that is on the right side of the backslash
                    species_name = re.search(r'(.*)/(.*)', org_name).group(2)
                    # Download url is constructed by using the website archetype for NCBI directories
                    # Data is organized by using a series of directories named as the species binomial
                    # For more information: https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/002/287/175/GCF_002287175.1_ASM228717v1/README.txt
                    #If FASTA genomic sequences are desired instead, change _protein.faa.gz to genomic.fna.gz
                    url_new = link + '/' + species_name + '_protein.faa.gz'
                    #If the protein file for a genome already exists, print a blank space
                    if os.path.exists(os.getcwd() + "/" + dom + "_protein_file" + "/" + species_name + "_protein.faa"):
                        print("")
                    else:
                        #Downloads faa file using the url constructed above
                        wget.download(url_new)
                        #Creates filename using genome accession number and ASM name
                        download_name = species_name + '_protein.faa.gz'
                        #Finds the location that the file has been downloaded to
                        source = os.path.abspath(download_name)
                        # moves file into domain specified folder
                        shutil.move(source, destination)
        #Closes opened assemble summary document
        assembly_summary.close()
        print("file_extraction--success")
    #Returns the location of the protein file >>> desired pathway/timestamp_domain/protein_file
    return destination

# Extracting all files and deleting the zip folder
# Input is the location of protein files
# No real output, the function extracts the compressed files
# Code was appropriated from: https://gist.github.com/kstreepy/a9800804c21367d5a8bde692318a18f5
def file_management(dest):
    dir_name = 'x'
    #References those files that are compressed
    extension = ".gz"
    file_name = ''
    os.chdir(dest)
    # Iterates through the protein file folder, destination where the .faa files were downloaded
    for item in os.listdir(dest):
        # check for ".gz" extension, does not apply to those files that were already extracted
        if item.endswith(extension):
            # get full path of files
            gz_name = os.path.abspath(item)
            # get file name for file within, GCF and ASM names
            file_name = (os.path.basename(gz_name)).rsplit('.', 1)[0]
            with gzip.open(gz_name, "rb") as f_in, open(file_name, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
            os.remove(gz_name)  # delete zipped file

    print("file_management-- success")


# DIAMOND IMPLEMENTATION: creates a new folder for DIAMOND outputs amd moves all outputs into the folder, creates a
# diamond-formatted library, completes DIAMOND search using a curated Uniprot library (reference.dmnd) as well as the
# protein files downloaded for each completely assembled genome.
# Input requires the location of unzipped protein files that was returned in file_extractor, as well as the reference document
# Output results in DIAMOND search results, which are lists of EC numbers found cited in each organisms protein sequence
# Function returns the location of the new folder, DIAMOND matches
#diamond_results_loc = diamond_impl(desired_location)
def diamond_impl(dest, name):
    os.chdir(dest)
    print(os.getcwd())
    # Creates another folder named DIAMOND matches to store DIAMOND output
    #new_destination = dest + '/DIAMOND_matches'
    # Issues with file management and moving to appropriate places, construct code after the distance matrix comes out right
    #Checks to see if the folder already exists, and if not creates a new folder
    # if os.path.exists(new_destination):
    #     print("hey i kicked myself out")
    #     #os.chdir(dest+'/DIAMOND_matches')
    # else:
    #     os.makedirs('DIAMOND_matches')
    # checks to see if the DIAMOND library is present in the new destination
    if os.path.exists('/home/anna/PycharmProjects/pythonProject/reference.dmnd'):
            #dest + '/reference.dmnd'):  # if there currently is no reference library (.dmnd), makedb creates a DIAMOND library
        print("Library Detected")
    # If not present, created a DIAMOND library by referencing the exact location where the Uniprot library is saved
    else:
        print("Creation of DIAMOND-formatted library...")
        makedb = ['diamond', 'makedb', '--in', '/home/anna/PycharmProjects/pythonProject/uniprot.fasta', '-d',
                  'reference']  # reference library full pathway
        subprocess.run(makedb)
        print("Library complete")
   # for item in domain_folder, if item is a faa file, complete diamond analysis
    for item in os.listdir(dest):
        # Checks for file extension
        if item.endswith('.faa'):
            #Finds path of file
            file_path = os.path.abspath(item)
            print(file_path)
            #finds the GCF/ASM name of the file by looking at the first part of the name before the .faa notation
            if name=="":
                file_name = (os.path.basename(file_path)).rsplit('.', 1)[0]
            else:
                file_name = name
            #New filename that ends with matches
            matches = file_name + "_matches.tsv"
            print(matches)
            # If genome has not already undergone DIAMOND search and is currently located in the correct folder
            # the subprocess function will run the diamond search
            new_destination = dest + "/" + matches
            if not path.isfile(matches):
                print("Processing ", file_name)
                #DIAMOND search using the full pathway of the protein files, max target sequence outputs only one best match with highest e-value
                blastp = ['diamond', 'blastp', '-d', 'reference.dmnd', '-q', file_path, '-o', matches,
                      '--max-target-seqs', '1', '--outfmt', '6']
                time.sleep(4)
                subprocess.run(blastp)
                # blastp = ['diamond', 'blastp', '-d', dest + '/reference.dmnd', '-q', file_path, '-o', matches,
                #           '--max-target-seqs', '1', '--outfmt', '6']
                #shutil.move(os.path.abspath(matches), new_destination)
                # else:
                #     #Moves DIAMOND outputs to designated folder
                #     if not re.match(os.path.abspath(matches), new_destination):
                #         shutil.move(os.path.abspath(matches), new_destination)
    print("diamond_impl--success")
    return dest