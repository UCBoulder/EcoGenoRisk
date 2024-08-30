from sklearn.cluster import KMeans
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

# np.random.seed(0)
# X = np.random.standard_normal((50,2)); #generates an array of shape 50 rows, 2 columns with a standard normal distribution 
# X[:25, 0] += 3 #modifies the first column, index 0, first 25 rows by adding 3 
# X[:25, 1] -= 4 #modifies the second column, index 1, first 25 rows by subtracting 4

#print(X)

full_data = pd.read_excel('/home/anna/Desktop/JD_Niche_OverLap (Git)/Niche_JD/Eco_V2/EnCen/Biome Analysis Results/River Sediment/Biome Synbio Top Match EC Comparison.xlsx')
def EC_parse(value):
    parts = (value.split('.'))
    return float('.'.join(parts[:2]))

full_data['EC Class and Subclass'] = full_data['EC Number'].apply(EC_parse)
# print(full_data.head(5))
# print(full_data['EC Class and Subclass'].dtype)


X = full_data.iloc[:, [1, 5]]

print(X)


# kmeans = KMeans(n_clusters=7, random_state=2, n_init=50).fit(X)
# print(kmeans.labels_)

# fig, ax = plt.subplots(1,1,figsize=(15,15))
f, ax = plt.subplots(figsize=(6.5, 6.5))
sns.despine(f, left=True, bottom=True)

sns.scatterplot(x="EC Class and Subclass", y="Percentage",
                hue="Synbio Presence/Absence",
                palette="ch:r=-.2,d=.3_r",
                sizes=(1, 8), linewidth=0,
                data=full_data, ax=ax)

plt.show()