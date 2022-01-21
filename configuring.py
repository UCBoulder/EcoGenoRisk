import wget
import requests
import re
import shutil
import os
import gzip
import subprocess
from subprocess import PIPE, Popen
from fake_useragent import UserAgent
import time

ua = UserAgent()
url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt'
path_folder = '/home/anna/PycharmProjects/pythonProject'
text_file = 'assembly_summary.txt'
header = {'User-Agent': str(ua.chrome)}

#Downloading/checking for assembly file
if os.path.exists(text_file):
     print("File exists")
else:
     temp_genome_list = requests.get(url, headers=header)
     temp_genome_list = requests.get(url)
     print(temp_genome_list.text)
     with open(text_file, 'w') as genome_list_out:
         genome_list_out.write(temp_genome_list.text)
     genome_list_out.close()

#looping through assembly summary file, finding the samples that have a complete genome, and downloading the FTP file
with open(text_file, 'r') as archaea_summary:
    temp= archaea_summary.readline()
    for line in archaea_summary:
        output = line.split("\t")
        if re.match(output[11],"Complete Genome"):
            link = output[19]
            org_name = link
            species_name = re.search(r'(.*)/(.*)', org_name).group(2)
            url_new = link+'/'+species_name+'_protein.faa.gz'
            #wget.download(url_new)
            download_name = species_name+'_protein.faa.gz'
            destination = os.path.abspath('archaea_protein_file')
            source = os.path.abspath(download_name)
            #shutil.move(source, destination)

#Extracting all files and deleting the zip folder
dir_name = 'x'
extension = ".gz"
file_name= ''
os.chdir(destination)
for item in os.listdir(destination):  # loop through items in dir
    if item.endswith(extension):  # check for ".gz" extension
        gz_name = os.path.abspath(item)  # get full path of files
        file_name = (os.path.basename(gz_name)).rsplit('.', 1)[0]  # get file name for file within
        with gzip.open(gz_name, "rb") as f_in, open(file_name, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
        os.remove(gz_name)  # delete zipped file


#DIAMOND IMPLEMENTATION
uniprot = open('/home/anna/PycharmProjects/pythonProject/uniprot.fasta')
reference = uniprot.read()
if os.path.abspath('home/anna/PycharmProjects/pythonProject/DIAMOND_matches') == 'EMPTY':
    directory = '/home/anna/PycharmProjects/pythonProject/DIAMOND_matches'
    os.makedirs(directory)


#SUBPROCESS INTERACTION
yeet = Popen(['diamond','version'])
print('stdout:', yeet.stdout)
makedb = ['diamond', 'makedb', '--', 'in', reference, '-d', 'reference']
lib=Popen(makedb, stdin=PIPE, stdout=PIPE)
#lib=subprocess.call(makedb)
lib.communicate()
counter = 0
for item in os.listdir(destination):
    if item.endswith('.faa'):
        print(counter)
        file_path = os.path.abspath(item)
        file_name = (os.path.basename(file_path)).rsplit('.', 1)[0]
        matches = file_name + "_matches"
        # run_diamond(file_name, reference, "blastp")
        time.sleep(10)
        blastp = ['diamond', 'blastp', '-d', 'reference', '-q', file_name, '-o', matches]
        time.sleep(10)
        #subprocess.run(str(blastp), capture_output=True, text=True)
        #search=Popen(blastp, stdin=PIPE,stdout=PIPE)
        #search=subprocess.call(blastp)
        #search.communicate()
        time.sleep(10)
        shutil.move(os.path.abspath(matches), directory)
        counter = counter + 1

archaea_summary.close()
uniprot.close()
print("That's all folks")
