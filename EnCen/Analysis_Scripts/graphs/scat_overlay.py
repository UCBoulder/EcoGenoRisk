#August 5nd, 2024
#This script shows the specific enzyme class and subclasses that are overlapped between the top match and the synbio org

import numpy as np 
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt
import re

activated_sludge = pd.read_excel('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/Paper Results/Biome Analysis Results/Activated Sludge/Biome Synbio Top Match EC Comparison.xlsx')
# print(activated_sludge.dtypes)

activated_sludge['EC subclass'] = activated_sludge['EC Number'].str.extract(r'^\d+\.(\d+)', expand = False)
activated_sludge['EC class'] = activated_sludge['EC Number'].str.extract(r'^(\d+)', expand = False) #converts the series into strings (previously objects), .extarct takes out specfic patterns, r denotes raw string, ^ asserts the position at the start of the string, \d matches any digit, () creates a group 
activated_sludge.to_excel('subclass.xlsx', index=False)
# print(activated_sludge.dtypes)

activated_sludge['EC subclass'] = activated_sludge['EC subclass'].astype(int)
activated_sludge['EC Class'] = activated_sludge['EC class'].astype(int)
# print(activated_sludge)

activated_sludge_tm = activated_sludge[['EC Class', 'EC subclass', 'Top Match Presence/Absence']]
print(activated_sludge_tm)
as_tm_present = activated_sludge_tm[activated_sludge_tm['Top Match Presence/Absence'] !=0.0]
activated_sludge_syn = activated_sludge[['EC Class', 'EC subclass', 'Synbio Presence/Absence' ]]
print(activated_sludge_syn)
as_syn_present = activated_sludge_syn[activated_sludge_syn['Synbio Presence/Absence'] !=0.0]

sns.scatterplot(data=as_syn_present, x = 'EC Class', y = 'EC subclass', color = 'red', label = 'Synbio')
sns.scatterplot(data=as_tm_present, x = 'EC Class', y = 'EC subclass', color = 'green', label = 'Top Match')
plt.subplots_adjust(left=0.15, right=0.85, top=0.85, bottom=0.15)  # Adjust margins
plt.title('Topmatch & synbio Enzyme Overlap')
plt.show()
