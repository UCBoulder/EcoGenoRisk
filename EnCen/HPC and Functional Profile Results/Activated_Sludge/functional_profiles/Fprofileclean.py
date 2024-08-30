##June 24th, 2024
##This script takes in the functional profile and turns it into a list of percentages within the biome, then merges those percentages to the input functional profle. Outputs percentages, EC's, EC names, and presence/abscence in a dataframe
import pandas as pd
import numpy as np
import requests
from fake_useragent import UserAgent
import time
from Bio import ExPASy
from Bio import SwissProt
from Bio.ExPASy import Enzyme
import pandas as pd 
import os
###______________________________________________________________________________________________________________
def percentages(biome_functional_profile):
    profile = pd.read_csv(biome_functional_profile, header=None)

    funcprof = profile.iloc[:,0].str.split(' ', expand=True) #this is the fancy bit of code that takes all the rows and first column (this is only one column in this dataframe to begin with) and splits each element by the space character. Expand = true puts each value in it's own column

    # print(funcprof.head())
    # print('total rows are: ', len(funcprof))
    # print('total columns are:', funcprof.shape[1])

    funcprof.columns = funcprof.iloc[0] #assigns the first row (iloc[0] to the columns of the dataframe)
    funcprof = funcprof[1:] #creates a new dataframe from the second row onward. 

    for col in funcprof.columns[1:]:
        funcprof[col] = pd.to_numeric(funcprof[col]) #converts all the dtype=object data to numbers. This is because pd.read_csv read the dataframe in as a string. 

    sum = funcprof.iloc[:, 1:].sum() #sums all rows, second column onward. 

    row = len(funcprof) -1 #just get's the amount of rows minus the summed row. 
    # print(row)
    funcprof.loc[len(funcprof)] = sum #adds the summed row onto the dataframe
    # print(funcprof)

    percentage = funcprof.iloc[1373].div(row) #divides the summed row by the amount of rows, i.e. finding the percentage of each enzyme. 
    print(percentage)
    final = pd.DataFrame(percentage) #turns numpy array into dataframe...why this is a numpy array, I don't know. 
    # final.reset_index(inplace=True)
    # print(final)
    final.columns = final.iloc[0] #turns the first row into the column names 
    final = final[1:] #new dataframe is second row onward. 
    # final = final.dropna()
    final.columns = ['Percentage'] #turns the column into a percentage 
    final = final.mul(100) #multiplies by 100 to get out of decimal form 
    # final = final.sort_values(by=['Percentage'], ascending=False) #sort the column percentage from high to low 
    final.reset_index(inplace=True)
    final.columns = ['EC Number', 'Percentage']
    
    return final
#At this point we have a dataframe of percentages and EC numbers______________________________________________________________________________________________________________________________________
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
###__________________________________________________________________________________________________________________________________________________________________________________________________
    handle = open(ec_library)
    records = Enzyme.parse(handle)
    ecnumbers = [record['ID'] for record in records]
    handle2 = open(ec_library)
    records = Enzyme.parse(handle2)
    ecnames = [record['DE'] for record in records]
    # print(type(ecnumbers)) #This is a list at this point in the code
    together = pd.DataFrame({'EC Number': ecnumbers, 'Name': ecnames})
    # together.to_csv('EC_name_&_num.csv', index=False)

    return together

def synbio_fun_profile(synbio_matrix):
    profile = pd.read_csv(synbio_matrix, header=None)

    funcprof = profile.iloc[:,0].str.split(' ', expand=True) #this is the fancy bit of code that takes all the rows and first column (this is only one column in this dataframe to begin with) and splits each element by the space character. Expand = true puts each value in it's own column
    # print(funcprof)
    # print(funcprof.head())
    # print('total rows are: ', len(funcprof))
    # print('total columns are:', funcprof.shape[1])

    # funcprof.columns = funcprof.iloc[0] #assigns the first row (iloc[0] to the columns of the dataframe)
    # funcprof = funcprof[1:] #creates a new dataframe from the second row onward. 
    # # print(funcprof)
    

    new_col = funcprof.iloc[0].tolist()
    new_col2 = funcprof.iloc[1].tolist()
    synbio_df = pd.DataFrame({'EC Number': new_col, 'Synbio Presence/Absence': new_col2})
    synbio_df = synbio_df.drop(0)
    # print(synbio_df)
    return synbio_df

def tm_fun_profile(synbio_matrix):
    profile = pd.read_csv(synbio_matrix, header=None)

    funcprof = profile.iloc[:,0].str.split(' ', expand=True) #this is the fancy bit of code that takes all the rows and first column (this is only one column in this dataframe to begin with) and splits each element by the space character. Expand = true puts each value in it's own column
    # print(funcprof)
    # print(funcprof.head())
    # print('total rows are: ', len(funcprof))
    # print('total columns are:', funcprof.shape[1])

    # funcprof.columns = funcprof.iloc[0] #assigns the first row (iloc[0] to the columns of the dataframe)
    # funcprof = funcprof[1:] #creates a new dataframe from the second row onward. 
    # # print(funcprof)
    

    new_col = funcprof.iloc[0].tolist()
    new_col2 = funcprof.iloc[1].tolist()
    synbio_df = pd.DataFrame({'EC Number': new_col, 'Top Match Presence/Absence': new_col2})
    synbio_df = synbio_df.drop(0)
    # print(synbio_df)
    return synbio_df

biome_functional_profile = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC Results/functional_profiles/Activated_Sludge_Metagenome_functional_profile'
final = percentages (biome_functional_profile)
together = EC_extract()
synbio_matrix = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC Results/functional_profiles/Synbio_functional_profile'
synbio = synbio_fun_profile(synbio_matrix)

top_match_bin = '/home/anna/Documents/UBA6164/UBA6164_functional_profile'
top_match = tm_fun_profile(top_match_bin)


#Now we have the dataframe of EC's and %'s, and the dataframe that has EC's and Names_______________________________________________________________________________________________________________
merged = pd.merge(final, together, how='inner', on='EC Number')
merged2 = pd.merge(merged, synbio, how='inner', on='EC Number')
merged3 = pd.merge(merged2, top_match, how='inner', on='EC Number')
merged4 = merged3.sort_values(by=['Percentage'], ascending=False) #sort the column percentage from high to low 
merged4.to_excel('Biome Synbio Top Match EC Comparison.xlsx', index=False)

