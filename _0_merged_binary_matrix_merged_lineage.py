import numpy as np
import pandas as pd
import os

##This documents purpose is to merge four datasets into one.
# (1) Archaea and Bacteria binary summary matrix must be combined into one summary matrix.
# (2) Archaea and Bacteria full lineage documents must be combined into one
# (3) Combined summary matrix and combined full lineage doc will be merged using the GCF names, under Name_of_Genome

##(1) Combining and formatting summary matrices
os.chdir('/home/anna/PycharmProjects/pythonProject')
#opens matrix summary files
bac_big_matrix = pd.read_csv('bacteria_combined1.csv', delimiter=" ", header=0, index_col=0)
arc_big_matrix = pd.read_csv('archaea_big_matrix.csv', delimiter=" ", header=0, index_col=0)
#Extracts only the name of genomes collected from NCBI
bac_index_for_editting = pd.DataFrame(bac_big_matrix.index, dtype='string')
arc_index_for_editting = pd.DataFrame(arc_big_matrix.index, dtype='string')
#Finds only the GCF Annotation by splitting data frames based on _ character
bac_editted=bac_index_for_editting.Name_of_Genome.str.split("_", expand= True)
arc_editted=arc_index_for_editting .Name_of_Genome.str.split("_", expand= True)
#Calls only on the column with GCF annotation
bac_big_matrix.index = 'GCF_'+ bac_editted.iloc[:,1]
arc_big_matrix.index = 'GCF_'+ arc_editted.iloc[:,1]
bac_big_matrix.index.name="Name_of_Genome"
arc_big_matrix.index.name="Name_of_Genome"
# Vertical concating of bacteria and archae EC binary matrix for hollistic analysis
bacteria_archaea = bac_big_matrix.append(arc_big_matrix)
# Outputs the total matrix dimensions
print(type(bacteria_archaea))
# Removes any previous annotations of the indeces
bacteria_archaea.index= bacteria_archaea.index.str.rstrip('_protein_matches')
print("Overall shape of archaea and bacteria summary matrix: ",np.shape(bacteria_archaea))
# Sets index column name as Name_of_Genome which is used for merging
#bacteria_archaea.to_csv('combined_binary_summary_matrix.csv' ,sep= '\t', header=True, index= True)
genomes_present = pd.DataFrame(bacteria_archaea.index, dtype='string')
#genomes_present.reset_index(inplace=True)
print(genomes_present)
genomes_present.to_csv('list_of_names.csv', sep='\t', header=True, index=False)

##(2) Combining and formatting of lineage documents
##LINEAGE DOCS ARE MANUALLY EDITTED ON LIBREOFFICE WITH THE FIND/REPLACE FUNCTION
# Initially downloaded documents are opened with tab and semi colon seperation to get the number of columns correctly
# The last updated time stamp is manually extracted from the website below
print('Last updated 2022-04-08 00:38 on GTDB: https://data.gtdb.ecogenomic.org/releases/latest/')
#Opens documents
bac_lineage = pd.read_csv('bacteria_lineage', delimiter='\t', dtype='string')
print(bac_lineage)
arc_lineage = pd.read_csv('archaea_lineage', delimiter='\t', dtype='string')
# Sets column names
bac_lineage.columns = ["Name_of_Genome", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
arc_lineage.columns = ["Name_of_Genome", "Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
# Vertically stacks the two lineage documents while maintaining only one column name
bacteria_archaea_lin = bac_lineage.append(arc_lineage)
bacteria_archaea_lin["Name_of_Genome"]= bacteria_archaea_lin["Name_of_Genome"].str.replace("+AF8-","_", regex = False)
print(bacteria_archaea_lin)
#bacteria_archaea_lin.reset_index(inplace=True)
# Prints out the total shape of the lneage key reference matrix
print("Lineage Key Matrix shape: ", np.shape(bacteria_archaea_lin))
bacteria_archaea_lin.to_csv('combined_lineage_doc.csv', sep='\t', header=True, index= True)
# Merges the combined binary matrix header with the lineage document
# This provides a full lineage document for the genomes we have collected from archaea and bacteria
# Note: this resultant matrix will be limited by the number of genomes we have collected
print('\n', genomes_present.columns,'\n',bacteria_archaea_lin.columns)

##(3) Merging of the two combined files based on GCF notation
#taxonomy_lineage = genomes_present.join(bacteria_archaea_lin.set_index(["Name_of_Genome"], verify_integrity=True), on= ["Name_of_Genome"])
taxonomy_lineage = pd.merge(genomes_present, bacteria_archaea_lin, on = ['Name_of_Genome'], how='inner')
#taxonomy_lineage = pd.merge(genomes_present, bacteria_archaea_lin, left_on = genomes_present.iloc[:,1], right_on=bacteria_archaea_lin.iloc[:,1])
#taxonomy_lineage = genomes_present.merge(bacteria_archaea_lin, how='inner', on = "Name_of_Genome")
#taxonomy_lineage = genomes_present.merge(bacteria_archaea_lin[list("Name_of_Genome")])
#taxonomy_lineage = genomes_present.set_index('Name_of_Genome').join(bacteria_archaea.set_index('Name_of_Genome'))
print("Total species identified: ", np.shape(taxonomy_lineage))
taxonomy_lineage.to_csv('taxonomy.csv' ,sep= '\t', header=True, index= True)
