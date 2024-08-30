import pandas as pd
import numpy as np 
from scipy import stats


# def genome_to_genome_diffcomp(domain_binary, comparison_domain_binary):
#     names_of_orgs = pd.DataFrame(domain_binary.index) 
#     diff = pd.DataFrame(abs(domain_binary.values - comparison_domain_binary.values)) 
#     row_sum = diff.sum(axis=1)
#     df1 = pd.DataFrame(row_sum)
#     difference_based_comparison = pd.concat([names_of_orgs, df1], axis=1)
#     difference_based_comparison.columns = ['Metagenomic Bin ID', 'Difference Score']
#     difference_based_comparison = difference_based_comparison.sort_values(by='Difference Score',
#                                                                           ignore_index=True).reset_index(drop=True)
#     difference_based_comparison.to_csv('Absolute_Difference_Comparison_Score.txt', header=True, index=True, sep='\t')




# domain_binary = pd.read_csv('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC Results/Activated_Sludge/functional_profiles/Activated_Sludge_Metagenome_functional_profile',
#                                   delimiter=" ", header=0)
# comparison_domain_binary = pd.read_csv('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC Results/Agricultural_Bulk_Soil/Agricultural_Bulk_Soil_Metagenome_functional_profile', delimiter=" ", header=0)

# genome_to_genome_diffcomp(domain_binary, comparison_domain_binary)

##Take the difference scores from two dataframes 
abs_as = pd.read_csv('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC Results/Activated_Sludge/synbio_inputs_and_outputs/Act_Sludge_Absolute_Difference_Comparison_Score.txt', delimiter='\t', header=0)
abs_agb = pd.read_csv('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC Results/Agricultural_Bulk_Soil/Ag_Bulk_soilAbsolute_Difference_Comparison_Score.txt', delimiter='\t', header=0)
abs_lake = pd.read_csv('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC Results/Lake_Sediment/Lake_Sediment_Absolute_Difference_Comparison_Score.txt', delimiter='\t', header=0 )
abs_river = pd.read_csv('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/HPC Results/River_Sediment/River_Sediment)Absolute_Difference_Comparison_Score.txt', delimiter='\t', header=0 )



# print(abs_as.head(5))
# as_list = abs_as['Difference Score'].tolist()
# abs_list = abs_agb['Difference Score'].tolist()

# Wilcox_result = stats.ranksums(as_list, abs_list)

# print("Statistic:", Wilcox_result.statistic) #The number 
# print("P-value:", Wilcox_result.pvalue)

lake_list = abs_lake['Difference Score'].tolist()
river_list = abs_river['Difference Score'].tolist()

Wilcox_result = stats.ranksums(lake_list, river_list)

print("Statistic:", Wilcox_result.statistic) #The number 
print("P-value:", Wilcox_result.pvalue)