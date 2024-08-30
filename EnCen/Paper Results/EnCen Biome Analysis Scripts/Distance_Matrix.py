###August 7th, 2024
##This script should take in the functional profiles and turn them into a di
import numpy as np 
import pandas as pd
from scipy.spatial.distance import pdist, squareform, cdist
import seaborn as sns
import matplotlib.pyplot as plt

def ec_sums(biome_functional_profile):
    profile = pd.read_csv(biome_functional_profile, header=None)

    funcprof = profile.iloc[:,0].str.split(' ', expand=True) #this is the fancy bit of code that takes all the rows and first column (this is only one column in this dataframe to begin with) and splits each element by the space character. Expand = true puts each value in it's own column

    # funcprof.to_excel('Funcprof1.xlsx', index=False)
    # print('total rows are: ', len(funcprof))
    # print('total columns are:', funcprof.shape[1])

    funcprof.columns = funcprof.iloc[0] #assigns the first row (iloc[0] to the columns of the dataframe)
    funcprof = funcprof[1:] #creates a new dataframe from the second row onward. 
    
    # funcprof.to_excel('Funcprof2.xlsx', index=False)

    for col in funcprof.columns[1:]:
        funcprof[col] = pd.to_numeric(funcprof[col]) #converts all the dtype=object data to numbers. This is because pd.read_csv read the dataframe in as a string. 

    # funcprof.to_excel('Funcprof3.xlsx', index=False)

    sum = funcprof.iloc[:, 1:].sum()

    return sum

def distance_matrix(as_totals, abs_totals):
    as_totals = as_totals.to_frame()
    abs_totals = abs_totals.to_frame()


    # print(as_totals)
    labels = ['Total']
    as_totals.columns = labels
    abs_totals.columns = labels
    as_totals.dropna(subset=['Total'])

    # print(as_totals)

    values1 = as_totals['Total'].values.reshape(-1, 1)
    values2 = abs_totals['Total'].values.reshape(-1, 1)

    distances = cdist(values1, values2, metric='euclidean')
    # distance_matrix = squareform(distances)
    result = pd.DataFrame(distances, index=as_totals.index, columns=abs_totals.index)

    return result


as_functional_profile = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC Results/Activated_Sludge/functional_profiles/Activated_Sludge_Metagenome_functional_profile'
abs_functional_profile = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC Results/Agricultural_Bulk_Soil/Agricultural_Bulk_Soil_Metagenome_functional_profile'
as_totals = ec_sums(as_functional_profile)
abs_totals = ec_sums(abs_functional_profile)

matrix = distance_matrix(as_totals, abs_totals)
print(matrix.index)




# sns.heatmap(matrix, cmap='coolwarm')
# plt.show()








