import pandas as pd
import wget
import requests
import re
import shutil
import os.path
import os
import protein_file
import subprocess
from fake_useragent import UserAgent

ua = UserAgent()
url = 'https://ftp.ncbi.nlm.nih.gov/genomes/refseq/archaea/assembly_summary.txt'
path_folder = '/home/anna/PycharmProjects/pythonProject'
text_file = 'assembly_summary.txt'
header = {'User-Agent': str(ua.chrome)}

if os.path.exists(text_file):
     print("File exists")
else:
     temp_genome_list = requests.get(url, headers=header)
     temp_genome_list = requests.get(url)
     print(temp_genome_list.text)
     with open(text_file, 'w') as genome_list_out:
         genome_list_out.write(temp_genome_list.text)
     genome_list_out.close()

processes = []
count = 0
total = 0

file_object = open(text_file, 'a')
with open(text_file ,'r') as archaea_summary:
    for line in archaea_summary:
        line1 = archaea_summary.readline()
        print(line1)
        if line1.find('Complete Genome') != -1:
            output = line1.split("\t")
            link = output[19]
            print(link)
            org_name = link
            species_name = re.search(r'(.*)/(.*)', org_name).group(2)
            url_new = link+'/'+species_name+'_protein.faa.gz'
            #print(line)
            #print(url_new)
            # wget.download(url_new)
            # file_name = species_name+'_protein.faa.gz'
            # destination = os.path.abspath('protein_file')
            # source = os.path.abspath(file_name)
            # shutil.move(source, destination+'/'+file_name)
            total+=1
        count+=1


print(count)
print(total)
archaea_summary.close()
#os.remove('assembly_summary.txt')
print("That's all folks")
