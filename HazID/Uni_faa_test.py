#####-----------------------------------------------------------###
#This script converts a 2 \t seperated TSV file into a FASTA file 
#Will then be plugged into _2_ line 200 
####------------------------------------------------------------###
import requests
import time
from fake_useragent import UserAgent
import os
import re
##_______________________________________________________________###
def tsv_to_fasta():
    reference_library = 'uniprot.tsv' 
    ua = UserAgent()
    header = {'User-Agent': str(ua.chrome)}
    uniprot_url = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cec%2Csequence&format=tsv&query=%28%28ec%3A*%29%29+AND+%28reviewed%3Atrue%29'
    time.sleep(4)
    uniprot = requests.get(uniprot_url, headers=header)
    if uniprot.status_code == 200:
        with open(reference_library, 'w+') as reference_library:
            reference_library.write(uniprot.text)
    print('Uniprot Reference File Has Been Created')
    
    input = os.path.abspath('uniprot.tsv')
    output = 'uniprot.fasta'
    with open(input, 'r') as input_file, open(output, 'w') as output_file:
        header=next(input_file)
        for line in input_file:
            carrot = f'>{line}'
            new_row = re.sub(r'(.*?\t.*?)\t', r'\1\n', carrot, 1)
            new_row_with_tab_replaced_by_question_mark = new_row.replace("\t", "?")
            # Be careful with the symbols used for replacement, as they have to align with the genomic summary string parsing
            output_row_with_space_replaced_with_ampersand = new_row_with_tab_replaced_by_question_mark.replace("; ", "&")
            #new_row2 = re.sub(r'(.*?\t.*?)\t', r'$', new_row, 0)
            #this is regex, for the record. 
            output_file.write(output_row_with_space_replaced_with_ampersand)
    print('FASTA has been created from TSV and is named', output)

tsv_to_fasta()















##File_Creation##
##This portion imports a TSV file from Uniprot using Uniprot's API
###Uniprot's API is generated with columns as "EC Number" and "Sequence"
####API is also generated with as Format = TSV and Compressed = No, using the API with streaming endpoint. 
#reference_library = '/projects/jodo9280/EcoDr/EcoDr/uniprot.tsv' 
#ua = UserAgent()
#header = {'User-Agent': str(ua.chrome)}
#uniprot_url = 'https://rest.uniprot.org/uniprotkb/stream?fields=accession%2Cec%2Csequence&format=tsv&query=%28%28ec%3A*%29%29+AND+%28reviewed%3Atrue%29'
#time.sleep(4)
#uniprot = requests.get(uniprot_url, headers=header)
#if uniprot.status_code == 200:
    #with open(reference_library, 'w+') as reference_library:
        #reference_library.write(uniprot.text)
#print('Uniprot Reference File Has Been Created')

##________________________________________________________________________###
##File Editing##
###This portion reads in the TSV file downloaded from NCBI, removes the first line, then iterates through each line to add a > and move the sequence down a row. 
#input = os.path.abspath('uniprot.tsv')
#output = '/projects/jodo9280/EcoDr/EcoDr/uniprot.fasta'
#with open(input, 'r') as input_file, open(output, 'w') as output_file:
    #header=next(input_file)
    #for line in input_file:
        #carrot = f'>{line}'
        #new_row = re.sub(r'(.*?\t.*?)\t', r'\1\n', carrot, 1)
        #output_file.write(new_row)

#print('FASTA has been created from TSV and is named', output)