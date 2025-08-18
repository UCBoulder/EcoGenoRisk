#August 6th, 2024
#This is a script that takes in a metagenomic functional profile and (hopefully) returns the bins that contain a specific enzyme. 


import pandas as pd 
import numpy as np


##Data cleaning and processing____________________________________________________________
profile = pd.read_csv('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC Results/Activated_Sludge/functional_profiles/Activated_Sludge_Metagenome_functional_profile', header=None)
funcprof = profile.iloc[:,0].str.split(' ', expand=True)
funcprof.columns = funcprof.iloc[0] #assigns the first row (iloc[0] to the columns of the dataframe)
funcprof = funcprof[1:] #creates a new dataframe from the second row onward. 
    
    # funcprof.to_excel('Funcprof2.xlsx', index=False)

for col in funcprof.columns[1:]:
    funcprof[col] = pd.to_numeric(funcprof[col])

# print(funcprof.head(5))

##Bin Search
value = 1
enzyme = '4.1.1.116'
bin_rows = funcprof[funcprof[enzyme] == value]
bin_number = bin_rows['Name_of_Genome']


bin_rows.set_index(bin_rows.columns[0], inplace=True)
print(bin_rows)

# one = bin_rows.loc['3300005657_9_matches.tsv'] - bin_rows.loc['3300012956_93_matches.tsv']
# two = bin_rows.loc['3300005657_9_matches.tsv'] - bin_rows.loc['3300055001_174_matches.tsv']
# three = bin_rows.loc['3300005657_9_matches.tsv'] - bin_rows.loc['3300056788_377_matches.tsv']
# four = bin_rows.loc['3300012956_93_matches.tsv'] - bin_rows.loc['3300055001_174_matches.tsv']
# five = bin_rows.loc['3300012956_93_matches.tsv'] - bin_rows.loc['3300056788_377_matches.tsv']
# six = bin_rows.loc['3300055001_174_matches.tsv'] - bin_rows.loc['3300056788_377_matches.tsv']
# # print(firstminussecond)

# sum1 = abs(pd.DataFrame(one).sum(axis=0))
# sum2 = abs(pd.DataFrame(two).sum(axis=0))
# sum3 = abs(pd.DataFrame(three).sum(axis=0))
# sum4 = abs(pd.DataFrame(four).sum(axis=0))
# sum5 = abs(pd.DataFrame(five).sum(axis=0))
# sum6 = abs(pd.DataFrame(six).sum(axis=0))



# # names = bin_rows['Name_of_Genome']
# # final = pd.DataFrame({'Bin': names, })
# results = [sum1, sum2, sum3, sum4, sum5, sum6]
# final_df = pd.DataFrame(results)
# # print(final_df)
