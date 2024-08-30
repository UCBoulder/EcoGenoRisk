###August 2nd, 2024: This script takes in, cleans, and merges the absolute difference scoring dataframes from HPC Results and creates a bar chart from the results. 
###The lowest score is the highest risk, reminder. 

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




activated_sludge = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC and Functional Profile Results/Activated_Sludge/synbio_inputs_and_outputs/Act_Sludge_Absolute_Difference_Comparison_Score.txt'
activated_sludge_df = abs_process(activated_sludge)
activated_sludge_df['Biome'] = 'Activated Sludge'

ag_bulk = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC and Functional Profile Results/Agricultural Bulk Soil/Bulk_Soil_Absolute_Difference_Comparison_Score.txt'
ag_bulk_df = abs_process(ag_bulk)
ag_bulk_df['Biome'] = 'Agricultural Bulk Soil'

cow_rumen = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC and Functional Profile Results/Cow Rumen/Cow_Rumen_Absolute_Difference_Comparison_Score.txt'
cow_rumen_df = abs_process(cow_rumen)
cow_rumen_df['Biome'] = 'Cow Rumen'

lake_sediment = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC and Functional Profile Results/Lake Sediment/Lake_Sediment_Absolute_Difference_Comparison_Score.txt'
lake_sediment_df = abs_process(lake_sediment)
lake_sediment_df['Biome'] = 'Lake Sediment'

river_sediment = '/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC and Functional Profile Results/River_Sediment/River_Sediment)Absolute_Difference_Comparison_Score.txt'
river_sediment_df =  abs_process(river_sediment)
river_sediment_df['Biome'] = 'River Sediment'



###Barplot___________________________________________________________________________________________________________________________________________________________________________________________
as5 = activated_sludge_df.head(5)
ab5 = ag_bulk_df.head(5)
cw5 = cow_rumen_df.head(5)
ls5 = lake_sediment_df.head(5)
rs5 = river_sediment_df.head(5)

combined1 = pd.concat([as5, ab5, cw5, rs5, ls5], ignore_index=True)
# combined2 = pd.concat([activated_sludge_df, ag_bulk_df, cow_rumen_df, river_sediment_df, lake_sediment_df], ignore_index=True)
# print(combined)
sns.barplot(data=combined1, x='Metagenome Bin ID', y='Difference Score', hue='Biome', palette='pastel', dodge=False)
plt.xticks(rotation=60)
plt.show()
# sns.barplot(data=combined2, x='Metagenome Bin ID', y='Difference Score', hue='Biome', palette='pastel', dodge=False)
# plt.xticks(range(len(combined['Metagenome Bin ID'].unique())), combined['Metagenome Bin ID'].unique())
# plt.xticks(rotation=60)
# plt.show()

###Swarmplot______________________________________________________________________________________________________________________________________________________________________________________________
sns.set_theme(style="whitegrid", palette="muted")
combined = pd.concat([activated_sludge_df, ag_bulk_df, cow_rumen_df, river_sediment_df, lake_sediment_df], ignore_index=True)

sns.swarmplot(data=combined, x='Difference Score', y='Biome', hue = 'Biome', s=2, legend = False)
plt.show()