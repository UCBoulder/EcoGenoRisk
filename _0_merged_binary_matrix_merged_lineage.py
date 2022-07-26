import numpy as np
import pandas as pd
import os as os

# This code's purpose is to merge four datasets into one.
# (1) Archaea and Bacteria binary summary matrix must be combined into one summary matrix.
# Binary summary matrix is the output from _3_genomic summary that produces an absence/presence matrix of enzyme
# commission (EC) numbers for every genome sampled
# (2) Archaea and Bacteria full lineage documents must be combined into one
# Lineage documents contain the full rank (Kingdom, Phylum, Class, Order, Family, Genus, Species) for all known Archaea
# and Bacteria organism (downloaded from Genome Taxonomy Database)
# (3) Combined summary matrix and combined full lineage doc will be merged using the assembly accession names for RefSeq
# assemblies (GCF names), under Name_of_Genome (NCBI Genome Assembly Model)
##====================================================================================================================##
# (1) Combining and formatting summary matrices
# Opens matrix summary files
bac_binary_path = '/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/Bacteria/bacteria_combined1.csv'
arc_binary_path = '/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/Archaea/archaea_big_matrix.csv'
bac_big_matrix = pd.read_csv(bac_binary_path, delimiter=' ', header = 0, index_col=0)
arc_big_matrix = pd.read_csv(arc_binary_path, delimiter=' ', header=0, index_col=0)
# Extracts only the name of genomes collected from National Center of Biotechnology Information (NCBI)
bac_index_for_editting = pd.DataFrame(bac_big_matrix.index, dtype='string')
arc_index_for_editting = pd.DataFrame(arc_big_matrix.index, dtype='string')

#Finds only the GCF Annotation by splitting data frames based on "_" character
bac_editted=bac_index_for_editting.Name_of_Genome.str.split('_', expand= True)
arc_editted=arc_index_for_editting .Name_of_Genome.str.split('_', expand= True)

#Calls only on the column with the GCF annotation
bac_big_matrix.index = 'GCF_'+ bac_editted.iloc[:,1]
arc_big_matrix.index = 'GCF_'+ arc_editted.iloc[:,1]
bac_big_matrix.index.name='Name_of_Genome'
arc_big_matrix.index.name='Name_of_Genome'

# Vertical concating of bacteria and archaea EC binary matrix to create an diverse overview of the two Kingdoms
bacteria_archaea = bac_big_matrix.append(arc_big_matrix)

# Prints the total matrix dimensions
print(type(bacteria_archaea))

# Removes any previous annotations of the indices
bacteria_archaea.index= bacteria_archaea.index.str.rstrip('_protein_matches')
print('Overall shape of archaea and bacteria summary matrix: ',np.shape(bacteria_archaea))

# Turn on if you would like to save the combined binary matrix summary for Archaea and Bacteria
bacteria_archaea.to_csv('combined_binary_summary_matrix.csv' ,sep= '\t', header=True, index= True)

# Saves as dataframe
genomes_present = pd.DataFrame(bacteria_archaea.index, dtype='string')
print(genomes_present)

#Turn on if you would like a list of all of the genomes collected from NCBI RefSeq
#genomes_present.to_csv('list_of_names.csv', sep='\t', header=True, index=False)
##====================================================================================================================##

# (2) Combining and formatting of lineage documents
# LINEAGE DOCS ARE MANUALLY EDITTED ON LIBREOFFICE WITH THE FIND/REPLACE FUNCTION
# Objects such as "-AF8-", "RS_", "GB_", "d_", "p_", "c_", "o_", "f_", "g_" were replaced manually
# Initially downloaded documents are opened with tab and semi colon seperation to get the correct number of columns
# The last updated time stamp is manually extracted from the website below
print('Last updated 2022-04-08 00:38 on GTDB: https://data.gtdb.ecogenomic.org/releases/latest/')

# Opens documents
bac_lineage = pd.read_csv('/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/Bacteria/bacteria_lineage', delimiter='\t', dtype='string')
print(bac_lineage)
arc_lineage = pd.read_csv('/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/Archaea/archaea_lineage', delimiter='\t', dtype='string')

# Sets column names
bac_lineage.columns = ['Name_of_Genome', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
arc_lineage.columns = ['Name_of_Genome', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

# Vertically stacks the two lineage documents while maintaining only one column name
bacteria_archaea_lin = bac_lineage.append(arc_lineage)
bacteria_archaea_lin['Name_of_Genome']= bacteria_archaea_lin['Name_of_Genome'].str.replace('+AF8-','_', regex = False)
print(bacteria_archaea_lin)

# Prints out the total shape of the lneage key reference matrix
print('Lineage Key Matrix shape: ', np.shape(bacteria_archaea_lin))

# Merges the combined binary matrix header with the lineage document
# This provides a full lineage document for the genomes collected for Archaea and Bacteria
# Note: this resultant matrix will be limited by the number of genomes
print('\n', genomes_present.columns,'\n',bacteria_archaea_lin.columns)

# (3) Merging of the two combined files based on GCF notation, keeps all of the entries based on the left input (genomes
# from NCBI), and replaces all of the blanks with NAs
taxonomy_lineage = pd.merge(genomes_present, bacteria_archaea_lin, on = ['Name_of_Genome'], how='left')
taxonomy_lineage.fillna('NA', inplace=True)
print('Total species identified: ', np.shape(taxonomy_lineage))
taxonomy_lineage.to_csv('taxonomy.tsv' ,sep= '\t', header=True, index= True)
print('Documents have been created for synbio analysis')
##===========================================Citations================================================================##
# NCBI Genome Assembly Model. https://www.ncbi.nlm.nih.gov/assembly/model/. Accessed 18 July 2022.