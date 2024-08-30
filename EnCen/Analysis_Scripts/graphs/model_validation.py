###August 15th, 2024
###This is for model validation, per https://www.sciencedirect.com/science/article/pii/S0048969720316405#s0085

import numpy as np 
import pandas as pd 
import seaborn as sns
import matplotlib.pyplot as plt

def abs_process(dataframe):
    ads = pd.read_csv(dataframe, sep='\t')
    # print(ads.head())
    ads2 = ads.drop(axis=1, columns = 'Unnamed: 0')
    ads2['Metagenome Bin ID'] = ads2['Metagenome Bin ID'].str.replace('_matches.tsv', '')
    # print(ads2.head())
    return ads2




Ag_Bulk = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/Paper Results/Model Validation/Bulk Soil R16 output/Bulk Soil R16/synbio_inputs_and_outputs/Absolute_Difference_Comparison_Score.txt'
Ag_Bulk_df = abs_process(Ag_Bulk)
Ag_Bulk_df['Biome'] = 'Agricultural Bulk Soil'

Rhizosphere = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/Paper Results/Model Validation/Soil Rhizosphere Outputs/OneDrive_2024-08-15(2)/Soil Rhizosphere/synbio_inputs_and_outputs/Absolute_Difference_Comparison_Score.txt'
Rhizosphere_df = abs_process(Rhizosphere)
Rhizosphere_df['Biome'] = 'Soil Rhizosphere'

Endosphere = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/Paper Results/Model Validation/Soil Endosphere Outputs/Absolute_Difference_Comparison_Score.txt'
Endosphere_df = abs_process(Endosphere)
Endosphere_df['Biome'] = 'Soil Endosphere'


print(type(Endosphere))

###Barplot___________________________________________________________________________________________________________________________________________________________________________________________
ab5 = Ag_Bulk_df.head(5)
rhizo5 = Rhizosphere_df.head(5)
endo5 = Endosphere_df.head(5)


combined = pd.concat([ab5, rhizo5, endo5], ignore_index=True)
combined = pd.concat([Ag_Bulk_df, Rhizosphere_df, Endosphere_df], ignore_index=True)
print(combined)
sns.barplot(data=combined, x='Metagenome Bin ID', y='Difference Score', hue='Biome', palette='pastel', dodge=False)
# plt.xticks(range(len(combined['Metagenome Bin ID'].unique())), combined['Metagenome Bin ID'].unique())
plt.xticks(rotation=60)
plt.show()