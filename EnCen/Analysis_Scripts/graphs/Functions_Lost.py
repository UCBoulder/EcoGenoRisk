##August 2nd, 2024
##This scirpt shows a bivariate blot of rare functions in the specified biome. 
##This script also shows a displot of the counts of each EC Plot. 

import numpy as np 
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt

activated_sludge= pd.read_excel('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/Paper Results/Biome Analysis Results/Cow Rumen/Lost_Functions.xlsx')
activated_sludge_all = pd.read_excel('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/Paper Results/Biome Analysis Results/Cow Rumen/Biome Synbio Top Match EC Comparison.xlsx')
 #remember lost functions are just those that the topmatch bin has and the invading organism does not have. 
# print(activated_sludge.head(10))

def EC_parse1(value):
    parts = (value.split('.'))
    return float('.'.join(parts[:2]))

def EC_parse2(value):
    parts = (value.split('.'))
    return float('.'.join(parts[:1]))

activated_sludge['EC Class and Subclass'] = activated_sludge['EC Number'].apply(EC_parse1)
activated_sludge2 = activated_sludge[['EC Class and Subclass', 'Percentage', 'Name']]
# activated_sludge2.to_excel('bivariate.xlsx', index=False)
# print(activated_sludge2)


##This is the Bivariate Plot__________________________________________________________________________________________________
# plot = sns.scatterplot(data=activated_sludge2, x= 'EC Class and Subclass', y='Percentage', color = 'blue', hue='Name',s=200) #Bivariate Plot 
# plt.tight_layout(rect=[0, 0, 0.75, 1])
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0, fontsize='large')
# plt.title('Percentages of Rare Functions Grouped by EC Class and Subclass')
# # plt.xticks([1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6, 6.5, 7.5])
# plot.xaxis.label.set_size(14)
# plot.yaxis.label.set_size(14)
# plt.show()

##This is the counts of each rare function_______________________________________________________________________________________
# activated_sludge['EC Class'] = activated_sludge['EC Number'].apply(EC_parse2)
# activated_sludge3 = activated_sludge[['EC Class', 'Percentage', 'Name']]
# print(activated_sludge3)
# sns.countplot(y='EC Class',  data=activated_sludge3)#hist, counts vs. EC
# plt.title('Counts of Rare Functions Grouped by EC Class')
# plt.show()


# activated_sludge_all = activated_sludge_all[activated_sludge_all['Percentage'] !=0.0]
# print(activated_sludge_all)

##This is visualizing the displaced functions in the entire funcitonal profile________________________________________________________
plt.figure(figsize=(12, 12))
sns.scatterplot(data= activated_sludge_all, x='Percentage', y='EC Number', color = 'grey', label = 'All Biome', legend=False)
plot = sns.scatterplot(data= activated_sludge, x='Percentage', y='EC Number', color='blue', hue = 'Name', palette='rocket', s=75, legend=False)
plt.subplots_adjust(left=0.15, right=0.85, top=0.85, bottom=0.15)  # Adjust margins
# plt.yticks(rotation=75)
y_values = activated_sludge['EC Number'].values
plt.yticks(ticks=[], rotation = 45)
plt.title('Metagenome Displaced Functions Overlayed on All Biome Functions')
plot.xaxis.label.set_size(14)
plot.yaxis.label.set_size(14)
plt.tick_params(left=False)
plt.show()