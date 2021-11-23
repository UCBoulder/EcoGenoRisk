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

print(header)
temp_genome_list = requests.get(url, headers=header)
with open(text_file,'w') as genome_list_out:
    genome_list_out.write(print(temp_genome_list))
genome_list_out.close()
processes = []
count =0
total=0

# if os.path.exists(text_file):
#     print("Assembly summary already exists")
# else:
#     wget.download(url, path_folder)

with open(text_file,'r') as archaea_summary:
    for line in archaea_summary:
        line1 = archaea_summary.readline()
        if "Complete Genome" in line1:
            output=line1.split("\t")
            link = output[19]
            print(link)
            #match= re.search("https://",line1)
            #print(match.group())
            #processes.append(line1.split()[-2])
            #print(processes[count])
            org_name = link
            species_name = re.search(r'(.*)/(.*)', org_name).group(2)
            # url_new = link[count]+'/'+species_name+'_protein.faa.gz'
            # print(url_new)
            # wget.download(url_new)
            # file_name = species_name+'_protein.faa.gz'
            # destination = os.path.abspath('protein_file')
            # source = os.path.abspath(file_name)
            # shutil.move(source, destination+'/'+file_name)
            count += 1
            total += 1
        else:
            total += 1
print(count)
print(total)
archaea_summary.close()
#os.remove('assembly_summary.txt')
print("That's all folks")

