import pandas as pd
import pubchempy as pcp
import numpy as np
import os as os
from itertools import chain


# For any PUGREST or PubChemPy issues, visit the following links:
# https://github.com/mcs07/PubChemPy/issues/33
# https://pubchempy.readthedocs.io/en/latest/api.html#search-functions
# https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest-tutorial

# Takes user input for chemical agent
# Connects to the PUGREST server and retrieves chemically similar compounds and their InChI-Key form
# Creates a summary dataframe called analogous
# In the case when the chemical agent does not have any similar/substructures/superstructures associated with its entry
# the program will return an "unavailable" tag to avoid PUGREST errors

# Ammonium: QGZKDVFQNNGYKY-UHFFFAOYSA-O
# Ammonia: QGZKDVFQNNGYKY-UHFFFAOYSA-N
# Carbon Dioxide: CURLTUGMZLYLDI-UHFFFAOYSA-N
# Butanol: LRHPLDYGYMQRHN-UHFFFAOYSA-N
def similarity_search(path, all_rxns_doc):
    os.chdir(path)
    print('Enter your inhibitors InChI-Key code: ')
    inhibitor = 'LRHPLDYGYMQRHN-UHFFFAOYSA-N'
    print('Enter number of similar structures desired: ')
    n = input()
    status = 1
    analogous = []
    try:
        inhibitor_info = pcp.get_compounds(inhibitor, 'inchikey', as_dataframe=True)
        inhibitor_cid = inhibitor_info.index[0].astype('str')
    except:
        print('Unable to find structure')
        status = 0
    try:
        print('Finding similar structures...')
        similar = pcp.get_compounds(inhibitor_cid, 'cid', searchtype='similarity', listkey_count=n)
        print(similar)
    except:
        print('Unable to find similar structures')
        status = 0
    try:
        print('Finding substructures...')
        substructure = pcp.get_compounds(inhibitor_cid, 'cid', searchtype='substructure', listkey_count=n)
        print(substructure)
    except:
        print('Unable to find substructures')
        status = 0
    try:
        print('Finding superstructures...')
        superstructure = pcp.get_compounds(inhibitor_cid, 'cid', searchtype='superstructure', listkey_count=n)
        print(superstructure)
    except:
        print('Unable to find superstructure')
        status = 0
    if status == 1:
        analogs = similar + substructure + superstructure
        inchikeys = []
        for compound in analogs:
            inchikeys.append(str(compound.inchikey))
        analogous = pd.DataFrame({'Similar_Structure': inchikeys})
        analogous.drop_duplicates()
        analogous.reset_index(drop=True)
        analogous.to_csv('similar_compounds.txt')
    return analogous, status


# Two inputs: dataset from BRENDA that contains the inhibitors and the organisms inhibited, and the InChI-Keys for all
# compounds given a common name or given a BioCyc CPD ID
# Function merges these datasets to translate the inhibitors into their InChI-Key names
# Output is a list of InChiKeys and EC numbers
# BRENDA dataset contains all organisms, downloaded on 2/15/2023
def inchikey_translation(inhibited_org, inchikey_compound_name):
    # Opens the BRENDA curated list of organisms and their inhibitors
    # List was created by searching for inhibitors for prokaryotes
    list_of_orgs_inhibited = pd.read_csv(inhibited_org, delimiter='\t')
    # Creates column names
    list_of_orgs_inhibited.columns = ['EC Number', 'Enzyme Name', 'Inhibitor Compound Name', '-', 'Species', '-', '-',
                                      'Primary accession Number']
    # Opens BioCyc curated list of compounds and their inchikeys
    comp_to_inchikey = pd.read_csv(inchikey_compound_name, delimiter='\t')
    comp_to_inchikey.columns = ['BioCyc Id', 'InChI-Key']
    # Converts the common names of compounds to all lower case (ex: CU 2+ --> cu 2+)
    list_of_orgs_inhibited['Inhibitor Compound Name'] = list_of_orgs_inhibited['Inhibitor Compound Name'].str.lower()
    comp_to_inchikey['BioCyc Id'] = comp_to_inchikey['BioCyc Id'].str.lower()
    # Merges two dataframes to create a list of cummulative list of compounds, their InChiKeys, and organisms inhibited
    inhibitors_translated = pd.merge(list_of_orgs_inhibited, comp_to_inchikey, how='inner',
                                     left_on='Inhibitor Compound Name',
                                     right_on='BioCyc Id')
    inhibitors_translated['InChI-Key'] = inhibitors_translated['InChI-Key'].str.replace('InChIKey=', '')
    # Saves the dataframe to a csv
    inhibitors_translated.to_csv('translated_inchikey.txt', header=True, index=True, sep='\t')
    # Subsets the data frame to select only the InChI-Key column and EC number
    # Overall dataframe is EC numbers that are inhibited along with their inhibitors InChI-Key notation
    ec_inchikeys_inhibited = inhibitors_translated[['EC Number', 'InChI-Key']]
    ec_inchikeys_inhibited.to_csv('ec_inchikey_inhibited.txt', sep='\t')
    return ec_inchikeys_inhibited


# Requires inputs from inchikey_translation(), similarity_search(), requires the EC BSM from NicheOverlap
# Using the combined EC BSM, EC numbers that are inhibited by the chemical agents and its analog are found
# List of all organisms that contain that EC number are also found
def vulnerable_pop(ec_inchikeys_inhibited, analogous):
    # Merges the EC number. inhibitor in InChiKey notation, and organism species name with the taxonomy dataset
    combined_binary = pd.read_csv(
        '/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/combined_binary_summary_matrix_2023_3_21.csv', header=0,
        index_col=0, sep='\t')
    # Creates a list of all EC numbers
    current_working_ec = pd.DataFrame(combined_binary.columns)
    # Creates column names
    current_working_ec.columns = ['Current Working ECs']
    # Merges analogous inhibitor list with list of EC numbers that are inhibited with their inhibitors
    # Overall result is EC number inhibited, Inhibitor Compound Name, Inhibitor InChI-Key name
    inhibited_by_sim_compounds_found = pd.merge(ec_inchikeys_inhibited, analogous, how='inner', left_on='InChI-Key',
                                                right_on='Similar_Structure')
    # Isolates only the EC number, InChiKey names
    inhibited_by_sim_compounds_abridged = inhibited_by_sim_compounds_found[['EC Number', 'InChI-Key']]
    # Saves the specified dataframe, isolates by column names
    inhibited_by_sim_compounds_abridged.to_csv('org_inhibited_in_GCF_notation.txt')
    # Creates a list of all unique EC numbers inhibited
    ec_numbers_inhibited = inhibited_by_sim_compounds_abridged['EC Number'].drop_duplicates()
    print('Total number of unique EC numbers inhibited is: ', ec_numbers_inhibited.shape)
    # Creates a list of EC numbers that are present in both the EC BSM and the inhibitors list in ec_numbers_inhibited
    # This addresses if new EC numbers are listed in the BRENDA dataset but not in the EC BSM
    ec_list_present = pd.merge(ec_numbers_inhibited, current_working_ec, how='inner', left_on='EC Number',
                               right_on='Current Working ECs')
    ec_list_present.to_csv('ec_inhibited.txt', sep='\t')
    # Creates summary dataframes
    for_particular_ec = pd.DataFrame(columns=['PoisInhibitor'])
    genomes = []
    # Loops through the EC numbers inhibited
    for ec in ec_list_present['EC Number']:
        print(ec)
        if ec == 'EC Number':
            print('')
        else:
            # If 1 is present in the EC column, isolate the rows with the present enzyme
            # Return a list of genomes that have the EC numebr present
            found_orgs = combined_binary[combined_binary.loc[:, ec] == 1].index.tolist()
            print('Number of genomes that contain the enzyme: ',
                  combined_binary[combined_binary.loc[:, ec] == 1].shape[0])
            if len(found_orgs) != 0:
                # Create a dataframe of all genomes that have the EC number inhibited
                new_df = pd.DataFrame(found_orgs, columns=[ec])
                # Vertically add the dataframe to the summary dataframe
                for_particular_ec = pd.concat([for_particular_ec, new_df], axis=1)
                # Add the genomes to a list
                genomes.append(found_orgs)
    # Opens the taxonomy spreadsheet
    taxonomy = pd.read_csv('/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/taxonomy_2023_3_21.tsv', header=0,
                           index_col=0, sep='\t')
    taxonomy.columns = ['Name_of_Genome', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    one_dim_genome = list(chain.from_iterable(genomes))
    df_genomes_inhibited = pd.DataFrame(one_dim_genome)
    df_genomes_inhibited.columns = ['Inhibited_Genomes']
    # Merges the GCF organisms with their taxonomic lineage
    labelled_inhibited_genomes = pd.merge(taxonomy, df_genomes_inhibited, left_on='Name_of_Genome',
                                          right_on='Inhibited_Genomes', how='inner')
    # Isolates the Species column and drops any duplicates, rows with NaN, empty rows, or rows with ,
    labelled_inhibited_genomes['Species'] = labelled_inhibited_genomes['Species'].drop_duplicates()
    labelled_inhibited_genomes.replace('', float('NaN'), inplace=True)
    labelled_inhibited_genomes.replace(',', float('NaN'), inplace=True)
    unique_labelled = labelled_inhibited_genomes['Species'].dropna(axis=0)
    unique_labelled.to_csv('labelled_inhibited_species.txt')
    return for_particular_ec, genomes


def scoring(inhibited_organisms_by_ec, genomes):
    # Converts list of all genomes in a numpy array
    as_array_genomes = np.array(genomes, dtype=object)
    # Finds the number of unique genomes
    unique_genomes = np.unique(as_array_genomes)
    print('Number of unique genomes is: ', unique_genomes.size)
    # Counts number of columns in the summary dataframe
    number_of_ec = inhibited_organisms_by_ec.shape[1] - 1
    print('Number of unique EC is: ', number_of_ec)
    inhibited_organisms_by_ec.to_csv('inhibited_organisms_by_ec.txt')
    return


path = '/home/anna/Desktop/EcoGenoRisk/HazID/PoisInbitor'
all_rxns_doc = '/home/anna/Desktop/EcoGenoRisk/HazID/CompetitorFind/All-reactions-of-MetaCyc.txt'
inhibited_org = '/home/anna/Desktop/EcoGenoRisk/HazID/PoisInbitor/prokaryote_inhibitor_list.csv'
inchikey_compound_name = '/home/anna/Desktop/EcoGenoRisk/HazID/CompetitorFind/InchiKeystoCompoundNames.txt'
[analogous, status] = similarity_search(path, all_rxns_doc)
if status == 1:
    ec_inchikeys_inhibited = inchikey_translation(inhibited_org, inchikey_compound_name)
    [for_particular_ec, genomes] = vulnerable_pop(ec_inchikeys_inhibited, analogous)
    scoring(for_particular_ec, genomes)
else:
    print('Try again with fewer number of similar structures or with a different InChI-Key')
