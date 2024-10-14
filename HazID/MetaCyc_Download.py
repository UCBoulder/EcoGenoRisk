#Imports

import pandas as pd
import numpy as np
import requests
import random

url =  'https://websvc.biocyc.org/st-get?id=biocyc14-86497-3916486221&format=tsv' #https://websvc.biocyc.org/st-get?id=[SMARTTABLE-ID]&format=[json|xml|tsv]
# Some notes onthis portion
# Metacyc smarttable MUST by set to public. Otherwise, the URL will break 
# Metacyc itself sucks to work with. Double check the URL has the column headings you need. Changes often do not save. 
params = {'random_param': random.randint(1, 1000)}
url = requests.get(url, params=params)
path = '/projects/jodo9280/EcoDr/EcoDr/All-reactions-of-MetaCyc.txt'
destination_file = '/projects/jodo9280/EcoDr/EcoDr/Competitor_Find/All-reactions-of-MetaCyc.txt'
if url.status_code == 200:
    with open(path, 'wb') as file:
        file.write(url.content)
    print("All Reactions of Metacyc Succesfully Downloaded")
else:
    print(f"Durin's Bane has Struck Once More. Status code: {url.status_code}")
MC = pd.read_csv('All-reactions-of-MetaCyc.txt', delimiter='\t')
MCEdit = MC[MC['Reaction'].str.contains('SUBSEQ,|\'end\'|(3)|is|beyond|the|end|of|the|sequence|(2).') == False]
#This gets rid of all the junk strings that come with the dataframe. 

MCEdit.to_csv('All-reactions-of-MetaCyc.txt', sep='\t', index=False, mode='w')
