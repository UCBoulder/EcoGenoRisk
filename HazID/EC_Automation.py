###_______________________________________________________________________________________________## 
#This is a test script for automating EC Number extraction. 
##This needs to take in a URL and extract just the ID from it, and put it into a CSV file so it can
##fit into the existing _3_ script
import requests
from fake_useragent import UserAgent
import time
from Bio import ExPASy
from Bio import SwissProt
from Bio.ExPASy import Enzyme
import os
import csv
###_______________________________________________________________________________________________##
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
    print(type(ec_file))
    print('EC List Has Been Created')
###_______________________________________________________________________________________________
    handle = open(ec_library)
    records = Enzyme.parse(handle)
    ecnumbers = [record["ID"] for record in records]
    print(type(ecnumbers)) #This is a list at this point in the code
    path = os.path.abspath(ec_library)
    with open(path, 'w+', newline='\n') as csv_file:
    #    csv_file = csv.writer(csv_file)
        for item in ecnumbers:
        #    csv_file.writerow([item])
            csv_file.write(item +'\n')
    print('EC list Has Been Created')
EC_extract()
