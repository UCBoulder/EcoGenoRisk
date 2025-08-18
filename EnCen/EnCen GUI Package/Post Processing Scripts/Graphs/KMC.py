#This script creates a KMeans clusters and visualizations 



from sklearn.cluster import KMeans
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import plotly.express as px
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
import matplotlib.pyplot as plt
import pandas as pd
from pandas.plotting import parallel_coordinates
from scipy.stats import hypergeom
from statsmodels.stats.multitest import multipletests


percentages = pd.read_excel('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/all_biomes_percentage.xlsx')

transposed = percentages.T
# print(transposed)
# dropped = percentages.drop(columns=['Biome'], errors='ignore')
transposed.columns = transposed.iloc[0]
transposed1 = transposed[1:]
# print(transposed1)

# num_clusters = 50
# kmeans_tests = [KMeans(n_clusters=i, init='random', n_init=10) for i in range(1, num_clusters)]
# score = [kmeans_tests[i].fit(transposed1).score(transposed1) for i in range(len(kmeans_tests))]

# # # Plot the curve
# plt.plot(range(1, num_clusters),score)
# plt.xlabel('Number of Clusters')
# plt.ylabel('Variance')
# # plt.show()

kmeans = KMeans(init='random', n_clusters=5, n_init=10, random_state=252)
transposed1['Cluster'] = kmeans.fit_predict(transposed1)


centroids = kmeans.cluster_centers_
print(centroids)

lables = kmeans.labels_
silhouette = silhouette_score(transposed1, lables)
# print(silhouette)
# print(transposed1)

fig = px.parallel_coordinates(transposed1)


colors = plt.cm.tab10.colors
plt.figure(figsize=(10, 6))
parallel_coordinates(transposed1, class_column='Cluster')
plt.show()


###This is graphing the kmeans outputs____________________________________________________________________________________________________________________________
c0 = transposed1[transposed1['Cluster'] == 0]
c1 = transposed1[transposed1['Cluster'] == 1]
c2 = transposed1[transposed1['Cluster'] == 2]
c3 = transposed1[transposed1['Cluster'] == 3]
c4 = transposed1[transposed1['Cluster'] == 4]
parallel_coordinates(transposed1, class_column='Cluster', color=colors)
plt.show()

plt.figure(figsize=(12, 8))
parallel_coordinates(c0, class_column='Cluster', color='black', alpha=0.25)
plt.ylim(0, 1.0)
plt.ylabel('EC Percentage')
plt.xticks(rotation=70)
plt.tight_layout()
plt.show()

plt.figure(figsize=(12, 8))
parallel_coordinates(c1, class_column='Cluster', color='black', alpha=0.25)
plt.ylim(0, 1.0)
plt.ylabel('EC Percentage')
plt.xticks(rotation=70)
plt.tight_layout()
plt.show()

plt.figure(figsize=(12, 8))
parallel_coordinates(c2, class_column='Cluster', color='black', alpha=0.25)
plt.ylim(0, 1.0)
plt.ylabel('EC Percentage')
plt.xticks(rotation=70)
plt.tight_layout()
plt.show()

plt.figure(figsize=(12, 8))
parallel_coordinates(c3, class_column='Cluster', color='black', alpha=0.25)
plt.ylim(0, 1.0)
plt.ylabel('EC Percentage')
plt.xticks(rotation=70)
plt.tight_layout()
plt.legend()
plt.show()

plt.figure(figsize=(12, 8))
parallel_coordinates(c4, class_column='Cluster', color='black', alpha=0.25)
plt.ylim(0, 1.0)
plt.ylabel('EC Percentage')
plt.xticks(rotation=70)
plt.tight_layout()
plt.show()


##Hypergeometric Analysis______________________________________________________________________________________________________________________

cluster = transposed1

cluster.sort_values("Cluster", axis=0, ascending=True, inplace=True)

reset = cluster.reset_index()

reset['EC'] = reset['index'].str.split('_').str[0]


reset2 = reset.drop(['index', 'Activated Sludge', 'Bulk Soil', 'Cow Rumen', 'Lake Sediment', 'River Sediment'], axis = 1)
reset2.columns.name  = None


reset2 = reset2.astype(int)


EC_per_cluster = reset2.groupby('Cluster')['EC'].value_counts()

print(EC_per_cluster)

total_EC_class = EC_per_cluster.groupby('EC').sum()
print(total_EC_class)

total_cluster_count = EC_per_cluster.groupby('Cluster').sum()
print(total_cluster_count)

# N=8370 #total population (Total EC's)
# K=2560 #Number of successes in the population (Total EC in the population)
# n=816 #Sample size (number of total ECs in the cluster)
# k=249 #Sample Successes (number of specific EC's in the cluster)

# p_value_1 = 1- stats.hypergeom.cdf(k-1, N, K, n)
# print(p_value_1)





#Cluster 0_________________________________________________________________________
N=8370 #total population (Total EC's) #Unchaged 
# K=2560 #Number of successes in the population (Total EC in the population) #changes
n=816 #Sample size (number of total ECs in the cluster) #Changes per cluster
# k=249 #Sample Successes (number of specific EC's in the cluster) #changes

K = [2560, 2287, 1806, 900, 347, 271, 98]
k = [249, 190, 188, 94, 44, 26, 25]

n2 = len(k)

def multi_p(N, n, K, k):
    p_values_final = []
    for i in range(n2):
        k2 = k[i]
        K2 = K[i]
        p_value = 1 - hypergeom.cdf(k2 -1, N, K2, n)
        # if p_value < 0.05: 
        p_values_final.append(p_value)
    return p_values_final

cluster_0_p_values = multi_p(N, n, K, k)
print('Cluster 0 p_values: ', cluster_0_p_values, '\n')


#Cluster 1_____________________________________________________________________________
n= 6720
k = [2187, 1845, 1472, 738, 259, 179, 40]


cluster_1_p_values = multi_p(N, n, K, k)
print('Cluster 1 p_values: ', cluster_1_p_values, '\n')

#Cluster 2_________________________________________________________________________________
n = 328
k = [68, 92, 86, 30, 18, 14, 20]


cluster_2_p_values = multi_p(N, n, K, k)
print('Cluster 2 p_values: ', cluster_2_p_values, '\n')

#Cluster 3___________________________________________________________________________________
n = 177
k = [16, 62, 31, 14, 15, 35, 4]


cluster_3_p_values = multi_p(N, n, K, k)
print('Cluster 3 p_values: ', cluster_3_p_values, '\n')

#Cluster 4_____________________________________________________________________________________
n = 228
k= [40, 98, 29, 24, 11, 17, 9]


cluster_4_p_values = multi_p(N, n, K, k)
print('Cluster 4 p_values: ', cluster_4_p_values, '\n')

#Bonferroni_Correction________________________________________________________________________________________
all_p_values = cluster_0_p_values + cluster_1_p_values + cluster_2_p_values + cluster_3_p_values +cluster_4_p_values
print(all_p_values)

b_rejected, b_corrected, _,_ = multipletests(all_p_values, alpha = 0.05, method='bonferroni')

print('rejected: ', b_rejected)
print('adjusted: ', b_corrected, '\n')

cluster_0_corrected = b_corrected[0:7]
print('Cluster 0 Correct: ', cluster_0_corrected, '\n')

cluster_1_corrected = b_corrected[7:14]
print('Cluster 1 Correct: ', cluster_1_corrected, '\n')


cluster_2_corrected = b_corrected[14:21]
print('Cluster 2 Correct: ', cluster_2_corrected, '\n')


cluster_3_corrected = b_corrected[21:28]
print('Cluster 3 Correct: ', cluster_3_corrected, '\n')


cluster_4_corrected = b_corrected[28:35]
print('Cluster 4 Correct: ', cluster_4_corrected, '\n')

print(len(all_p_values))

#Histograms________________________________________________________________________________________________________________________
#Cluster 0__________________
k0 = [249, 190, 188, 94, 44, 26, 25] 
total = 816

result0 = [ x / total for x in k0]
x = [1, 2, 3, 4, 5, 6, 7]
print(result0)


sorted_indices = np.argsort(result0)[::-1]
sorted_result0 = np.array(result0)[sorted_indices]
sorted_x = np.array(x)[sorted_indices]

# Create a grayscale color palette and reverse it
colors = sns.color_palette("Greys", len(sorted_result0))[::-1]

# Map the sorted values to the color palette
color_mapping = {sorted_result0[i]: colors[i] for i in range(len(sorted_result0))}
bar_colors = [color_mapping[val] for val in result0]

# Plot with original data and mapped colors
sns.barplot(x=x, y=result0, palette=bar_colors, errorbar='sd')
# plt.xlabel('Data Points')
# plt.ylabel('Proportion')
# plt.title('Cluster 0 Proportions')
plt.ylim(0, 0.5)
# plt.show()
print(np.var(result0))

# #Cluster 1______________________________________
k0 = [2187, 1845, 1472, 738, 259, 179, 40]
total = 6720

result0 = [ x / total for x in k0]
x = [1, 2, 3, 4, 5, 6, 7]
print(result0)


sorted_indices = np.argsort(result0)[::-1]
sorted_result0 = np.array(result0)[sorted_indices]
sorted_x = np.array(x)[sorted_indices]

# Create a grayscale color palette and reverse it
colors = sns.color_palette("Greys", len(sorted_result0))[::-1]

# Map the sorted values to the color palette
color_mapping = {sorted_result0[i]: colors[i] for i in range(len(sorted_result0))}
bar_colors = [color_mapping[val] for val in result0]

# Plot with original data and mapped colors
sns.barplot(x=x, y=result0, palette=bar_colors)
# plt.xlabel('Data Points')
# plt.ylabel('Proportion')
# plt.title('Cluster 0 Proportions')
plt.ylim(0, 0.5)
plt.show()
print(np.var(result0))



# #Cluster 2_____________________________________________
k0 = [68, 92, 86, 30, 18, 14, 20]
total = 328

result0 = [ x / total for x in k0]
x = [1, 2, 3, 4, 5, 6, 7]
print(result0)


sorted_indices = np.argsort(result0)[::-1]
sorted_result0 = np.array(result0)[sorted_indices]
sorted_x = np.array(x)[sorted_indices]

# Create a grayscale color palette and reverse it
colors = sns.color_palette("Greys", len(sorted_result0))[::-1]

# Map the sorted values to the color palette
color_mapping = {sorted_result0[i]: colors[i] for i in range(len(sorted_result0))}
bar_colors = [color_mapping[val] for val in result0]

# Plot with original data and mapped colors
sns.barplot(x=x, y=result0, palette=bar_colors)
# plt.xlabel('Data Points')
# plt.ylabel('Proportion')
# plt.title('Cluster 0 Proportions')
plt.ylim(0, 0.5)
plt.show()
print(np.var(result0))




# #Cluster 3_______________________________________________
k0 = [16, 62, 31, 14, 15, 35, 4]
total = 177

result0 = [ x / total for x in k0]
x = [1, 2, 3, 4, 5, 6, 7]
print(result0)

sorted_indices = np.argsort(result0)[::-1]
sorted_result0 = np.array(result0)[sorted_indices]
sorted_x = np.array(x)[sorted_indices]

# Create a grayscale color palette and reverse it
colors = sns.color_palette("Greys", len(sorted_result0))[::-1]

# Map the sorted values to the color palette
color_mapping = {sorted_result0[i]: colors[i] for i in range(len(sorted_result0))}
bar_colors = [color_mapping[val] for val in result0]

# Plot with original data and mapped colors
sns.barplot(x=x, y=result0, palette=bar_colors)
# plt.xlabel('Data Points')
# plt.ylabel('Proportion')
# plt.title('Cluster 0 Proportions')
plt.ylim(0, 0.5)
plt.show()
print(np.var(result0))



# #Cluster 4________________________________________________
k0 = [40, 98, 29, 24, 11, 17, 9]
total = 228

result0 = [ x / total for x in k0]
x = [1, 2, 3, 4, 5, 6, 7]
print(result0)

sorted_indices = np.argsort(result0)[::-1]
sorted_result0 = np.array(result0)[sorted_indices]
sorted_x = np.array(x)[sorted_indices]

# Create a grayscale color palette and reverse it
colors = sns.color_palette("Greys", len(sorted_result0))[::-1]

# Map the sorted values to the color palette
color_mapping = {sorted_result0[i]: colors[i] for i in range(len(sorted_result0))}
bar_colors = [color_mapping[val] for val in result0]

# Plot with original data and mapped colors
sns.barplot(x=x, y=result0, palette=bar_colors)
# plt.xlabel('Data Points')
# plt.ylabel('Proportion')
# plt.title('Cluster 0 Proportions')
plt.ylim(0, 0.5)
plt.show()
print(np.var(result0))
