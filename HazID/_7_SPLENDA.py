###______________________________________________________________________________________________________________________________________________###
#John Docter, March 11th, 2024. 
##This script eats BRENDA's JSON download and converts it into a usuable format for Pois Inhibitor

import numpy as np
import json
import pandas as pd
from datetime import datetime
import os as os

path = '/projects/jodo9280/EcoDr/EcoDr/Poisinhibitor'
with open('/projects/jodo9280/EcoDr/EcoDr/Poisinhibitor/brenda_2023_1.json', 'r') as file:
        data1 = json.load(file)

os.chdir(path) 
date = datetime.now()
date2 = date.strftime('%Y-%m-%d')
filename = 'BRENDA_Inhibitor_list: Updated {}.csv'.format(date2)
print('Beep Bop Boop.....Creating BRENDA Inhibitor List')

##Attempt 1__________________________________________________________________________________________________________________________________________________________________________
def SPLENDA(): 
    with open(filename, 'w') as f:
        for key, value in data1["data"].items():
        # Access the "cofactor" for each key
            print(key, file = f)
            # inhibitor_list.append(inhibitor)
            inhibitor = value.get("inhibitor")
            if inhibitor:
                for compound in inhibitor:
                    # print("Inhibitor:", compound)
                    value = compound.get("value")
                    if value:
                        print("Value:", value, file =f)
                    else:
                        print("No inhibitor value found")
    
    BRENDA1 = pd.read_csv(filename, sep = '\t', header=None)
    # print(BRENDA1)
    BRENDA1.columns = ['EC Number']
    BRENDA1['Inhibitor'] = ''
    mask = BRENDA1['EC Number'].str.contains('Value')
    # print(mask)
    # shifted = BRENDA1['EC'].shift(1)

    BRENDA1.loc[mask, 'Inhibitor'] = BRENDA1['EC Number'].shift(0)
    BRENDA1.loc[mask, 'EC Number'] = None
    BRENDA1['EC Number'] = BRENDA1['EC Number'].str.replace('spontaneous', '')
    BRENDA1.drop(0, inplace=True)
    BRENDA1['Inhibitor'] = BRENDA1['Inhibitor'].shift(-1)
    BRENDA1['EC Number']=BRENDA1['EC Number'].ffill()
    # (method='ffill', inplace=True)
    BRENDA1.replace(r'^\s*$', pd.NA, regex=True, inplace=True)
    BRENDA1.dropna(subset = ['Inhibitor'], inplace=True)
    BRENDA1['Inhibitor'] = BRENDA1['Inhibitor'].str.replace('Value: ','')
    BRENDA1.to_csv(filename, sep='\t')
    print('Bloop Inhibitor List...Complete')

SPLENDA()
