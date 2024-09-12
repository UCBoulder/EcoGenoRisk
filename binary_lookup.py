import pandas as pd 

synbio_binary = '/home/anna/Desktop/Synbio_functional_profile'


synbio_binary = pd.read_csv(synbio_binary, delimiter=" ", header=0)
# print(type(synbio_binary))
# synbio_binary = synbio_binary.set_index('Name_of_MetaGenome_Bin')
# print(synbio_binary)
transpose = synbio_binary.transpose()
# print(transpose) #fuckin bingo 
# print(transpose.index)
value = transpose.loc['2.3.1.1']
print(value)