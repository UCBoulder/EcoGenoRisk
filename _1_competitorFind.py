import numpy
import pandas as pd
import numpy as np
import os as os
from scoring import scoring, finding_distance
#nonsense add

##====================================================================================================================##

# Inputs are any data frame that needs transforming
# Converts any dataframe into a single column dataframe
# Outputs the single column df
def to_one_column(df):
    all_values = []
    for column in df:
        # Converts every line to a list and adds it to itself
        pathway_list = df[column].tolist()
        all_values += pathway_list
    one_col = pd.DataFrame(all_values, columns=['InChI-Key'])
    # Removes all None, blank spaces with NaN
    one_col['InChI-Key'].replace('None', np.nan)
    one_col['InChI-Key'].replace('', np.nan)
    # Drops all NaN types
    one_col.dropna()
    # Converts all lines into string and removes white space when present
    one_col['InChI-Key'].astype(str).str.strip()
    # Removes any duplicates
    one_col['InChI-Key'].drop_duplicates()
    # Resets the index in single column
    one_col.reset_index(drop=True)
    print(one_col)
    return one_col


##====================================================================================================================##

# Filters out universal InChI-Keys such as water, NADP, NADPH, Proton
# Input is a merged function from substrate_changes_synbio_v_chassis() or substrate_changes_modified_pathway()
# Uses document InchiKeystoCompoundNames.txt that was curated on MetaCyc using custom SmartTables. Document contains two
# columns: compound name, and Inchikey id for corresponding compound.

def relevant_compounds(df):
    # Finds the indexes where InChiKey for water are present
    index_water = df[df['InChI-Key'] == 'InChIKey=XLYOFNOQVPJJNP-UHFFFAOYSA-N'].index
    # Finds the indexes where InChiKey for NADPH are present
    index_NADPH = df[df['InChI-Key'] == 'InChIKey=ACFIXJIJDZMPPO-NNYOXOHSSA-J'].index
    # Finds the indexes where InChiKey for NADP are present
    index_NADP = df[df['InChI-Key'] == 'InChIKey=XJLXINKUBYWONI-NNYOXOHSSA-K'].index
    # Finds the indexes where InChiKey for Proton are present
    index_proton = df[df['InChI-Key'] == 'InChIKey=GPRLSGONYQIRFK-UHFFFAOYSA-N'].index
    # Removes row if InChIKEy for water is present
    df.drop(index_water, inplace=True)
    # Removes row if InChIKey for NADPH is present
    df.drop(index_NADPH, inplace=True)
    # Removes row if InChIKey for NADP is present
    df.drop(index_NADP, inplace=True)
    # Removes row if InChIKey for Proton is present
    df.drop(index_proton, inplace=True)
    # Resets index
    df.reset_index(drop=True)
    # Returns list of relevant InChiKeys
    return df


##====================================================================================================================##

# Translates InChI-Keys into conventional compound names
# Returns dataframe of Inchikeys and its translation
def inchikey_to_conventional_names(df):
    key = pd.read_csv('InchiKeystoCompoundNames.txt', delimiter='\t', header=0, index_col=None)
    translated_inchikeys = pd.merge(df, key, on='InChI-Key', how='left')
    return translated_inchikeys


##====================================================================================================================##

# Inputs for this function are the pathway where synbio analysis documents were saved, and the name of the synbio organism
# Outputs include: a data frame which contains the EC numbers that are variable between the two genomes and which
# organism the EC numbers occur, the top match's genome name and EC BSM
# This function tracks changes in the EC BSM between the synbio organism and the top match found in the function
# genome_to_genome_diffcomp() from the _5_synbio_distance script
def ec_locator(sb_name, pathway):
    # Reads the document that scores all genomes based on EC expression
    find_top_match = pd.read_csv(pathway + '/Chimera1_Difference_Based_Comparison_Score.txt', delimiter='\t',
                                 header=0,
                                 index_col=0)
    # Finds the top match row from the document. This organism is functionally most similar to the synbio organism
    top_match_row = find_top_match.head(1).reset_index(drop=True)
    # Extracts the top match genome name
    #top_match_name = top_match_row.iloc[0, 0]
    top_match_name = 'GCF_001953035.1'
    print('Top Match Found is: ', top_match_name)
    # Reads document which has EC BSM for Bacteria, Archaea, and the synbio organism
    ec_bsm = pd.read_csv(pathway + '/complete_binary_matrix.txt', delimiter='\t', header=0, index_col=0)
    # Finds the EC BSM for top match organism and transforms it into a data frame
    top_match_bsm = pd.DataFrame(ec_bsm.loc[top_match_name], columns=[top_match_name])
    # Sets index name of dataframe to EC Number
    top_match_bsm.index.name = 'EC-Number'
    # Finds the EC BSM for synbio organism and transforms it into a data frame
    synbio_bsm = pd.DataFrame(ec_bsm.loc[sb_name], columns=[sb_name])
    # Sets index name of dataframe to EC Number
    synbio_bsm.index.name = 'EC-Number'
    # Merges two dataframes based on the EC Number rows based on outer, creating a large combined dataframe with EC BSM
    # for synbio EC BSM and top match EC BSM
    different_ECs = pd.merge(top_match_bsm, synbio_bsm, on='EC-Number', how='outer')
    # Finds occurences where the synbio and top match EC BSMs don't have the same binary values and extracts EC numbers
    different_ECs = different_ECs[different_ECs[sb_name] != different_ECs[top_match_name]]
    # Turn on to save the list of EC numbers that are different between the two genomes and where the EC number is present
    different_ECs.to_csv(sb_name + '_different_EC_profile.txt', header=True, index=True, sep='\t')
    top_match_bsm = top_match_bsm[(top_match_bsm.loc[:, top_match_name] != 0)]
    synbio_bsm = synbio_bsm[(synbio_bsm.loc[:, sb_name] != 0)]
    print('Locating Variable EC numbers Complete')
    return different_ECs, top_match_bsm, synbio_bsm, top_match_name


##====================================================================================================================##

# Inputs for this function include: location where the CompetitorFind framework is saved,
# top matches original EC BSM, synbio EC BSM
# Output for this function is dataframe which contain InChIKeys of substrates that are shared between the synbio
# organism and the top match
# Method is a "catch all", no matter what kind of changes occur in the EC BSM
def substrate_changes_synbio_v_chassis(sb_name, pathway, chassis_bsm, synbio_bsm):
    # Opens MetaCyc list of all reactions
    metacyc_all_rxns = pd.read_csv(pathway + '/All-reactions-of-MetaCyc.txt',
                                   delimiter='\t', header=0, index_col=0)
    # Converts into database
    metacyc_all_rxns = pd.DataFrame(metacyc_all_rxns)
    print('Old column names', metacyc_all_rxns.columns)
    metacyc_all_rxns.columns = ['EC-Number', 'Substrates', 'Substrates InChI-Key', 'Reactants',
                                'Reactants InChI-Key', 'Products', 'Products InChI-Key']
    # String processing of dataframe by removing 'EC-' from the start of the EC number
    print('Here is the reactants column', metacyc_all_rxns['Reactants InChI-Key'])
    metacyc_all_rxns['EC-Number'] = metacyc_all_rxns['EC-Number'].str.replace('EC-', '', regex=False)
    # Merges based on EC number to create a list of reactions/substrates occurring in top match
    top_match_rxns = pd.merge(chassis_bsm, metacyc_all_rxns, on='EC-Number', how='inner')
    print(top_match_rxns.columns)
    print(top_match_rxns['Reactants InChI-Key'])
    # Turn on to save the list of substrates found in top match
    top_match_rxns.to_csv('topmatch_' + sb_name + '_chassis_all_rxns.txt', header=True, index=True, sep='\t')
    # Isolates the InChI Key column and splits the column based on // to isolate all substrates for top match
    chassis_InChI_Key = pd.DataFrame(top_match_rxns['Reactants InChI-Key'].astype(str).str.split('//', expand=True))
    chassis_one_col = to_one_column(chassis_InChI_Key)
    chassis_one_col.to_csv(sb_name + '_chassis_all_rxns.txt', header=True, index=True, sep='\t')
    # Merges based on EC number to create a list of reactions/substrates occurring in synbio
    synbio_rxns = pd.merge(synbio_bsm, metacyc_all_rxns, on='EC-Number', how='inner')
    # Merges the synbio binary summary matrix with the metacyc list to find the InChI-Key lists
    # Splits InChI-Keys based on the '//' separator
    synbio_InChI_Key = pd.DataFrame(synbio_rxns['Reactants InChI-Key'].astype(str).str.split('//', expand=True))
    synbio_one_col = to_one_column(synbio_InChI_Key)
    # Turn on to save the list of substrates found in synbio
    # synbio_one_col.to_csv(sb_name+'_synbio_all_rxns.txt', header=True, index= True, sep='\t')
    # Isolates the InChI Key column and splits the column based on // to isolate all substrates for synbio
    # Creates an array of InChI Keys of substrates that can be found in both organisms, saves the array as index
    shared_InChI_Key = pd.merge(synbio_one_col, chassis_one_col, on='InChI-Key', how='inner').reset_index(drop=True)
    print('I have merged')
    # Returns the InChI Key names by referencing the index
    # Converts array into a single column list
    # Removes any spaces
    shared_InChI_Key['InChI-Key'] = shared_InChI_Key['InChI-Key'].str.strip()
    # Finds unique InChI Keys in the list
    unique_chassis_InChI_Key = shared_InChI_Key['InChI-Key'].drop_duplicates()
    unique_chassis_InChI_Key.reset_index(drop=True)
    # Saves list of InChI Keys
    unique_chassis_InChI_Key = pd.DataFrame(unique_chassis_InChI_Key, columns=['InChI-Key'])
    # Removes the common InChI-Keys such as proton, ATP, and saves the list
    chassis_inchi_keys_translated = relevant_compounds(unique_chassis_InChI_Key)
    chassis_inchi_keys_translated.to_csv(sb_name + '_synbiovschassis_inchikey.txt', header=True, index=True, sep='\t')
    # unique_chassis_InChI_Key.to_csv(sb_name+'_synbiovschassis_inchikey.txt', header=True, index= True, sep='\t')
    print('Top Match vs. Synbio InChI Key Substrates Analysis Is Complete')
    return unique_chassis_InChI_Key


##====================================================================================================================##

# Inputs for this function include: dataframe which lists all the different EC numbers found in ec_locator(),
# the location where the CompetitorFind framework is saved
# Output for this function is dataframe which contain InChIKeys of substrates that are shared in the modified pathways
# of the synbio and top match, output for one is input for another = mutualism
# This method is not effective for unilateral changes (e.i. synbio is created by turning on EC numbers from top match)
# Method is preferred if synbio and top match both have EC numbers turned on and off

def mutualism1_modified_pathway(sb_name, pathway, chassis_bsm, synbio_bsm):
    # Opens MetaCyc list of all reactions
    metacyc_all_rxns = pd.read_csv(pathway + '/All-reactions-of-MetaCyc.txt',
                                   delimiter='\t', header=0, index_col=0)
    # Converts into database
    metacyc_all_rxns = pd.DataFrame(metacyc_all_rxns)
    print('Old column names', metacyc_all_rxns.columns)
    metacyc_all_rxns.columns = ['EC-Number', 'Substrates', 'Substrates InChI-Key', 'Reactants',
                                'Reactants InChI-Key', 'Products', 'Products InChI-Key']
    # String processing of dataframe by removing 'EC-' from the start of the EC number
    print('Here is the reactants column', metacyc_all_rxns['Reactants InChI-Key'])
    metacyc_all_rxns['EC-Number'] = metacyc_all_rxns['EC-Number'].str.replace('EC-', '', regex=False)
    # Merges based on EC number to create a list of reactions/substrates occurring in top match
    top_match_rxns = pd.merge(chassis_bsm, metacyc_all_rxns, on='EC-Number', how='inner')
    print(top_match_rxns.columns)
    print(top_match_rxns['Reactants InChI-Key'])
    # Turn on to save the list of substrates found in top match
    top_match_rxns.to_csv('topmatch_' + sb_name + '_chassis_all_rxns.txt', header=True, index=True, sep='\t')
    # Isolates the InChI Key column and splits the column based on // to isolate all substrates for top match
    chassis_InChI_Key = pd.DataFrame(top_match_rxns['Products InChI-Key'].astype(str).str.split('//', expand=True))
    chassis_one_col = to_one_column(chassis_InChI_Key)
    chassis_one_col.to_csv(sb_name + '_chassis_all_rxns.txt', header=True, index=True, sep='\t')
    # Merges based on EC number to create a list of reactions/substrates occurring in synbio
    synbio_rxns = pd.merge(synbio_bsm, metacyc_all_rxns, on='EC-Number', how='inner')
    # Merges the synbio binary summary matrix with the metacyc list to find the InChI-Key lists
    # Splits InChI-Keys based on the '//' separator
    synbio_InChI_Key = pd.DataFrame(synbio_rxns['Reactants InChI-Key'].astype(str).str.split('//', expand=True))
    synbio_one_col = to_one_column(synbio_InChI_Key)
    # Turn on to save the list of substrates found in synbio
    # synbio_one_col.to_csv(sb_name+'_synbio_all_rxns.txt', header=True, index= True, sep='\t')
    # Isolates the InChI Key column and splits the column based on // to isolate all substrates for synbio
    # Creates an array of InChI Keys of substrates that can be found in both organisms, saves the array as index
    shared_InChI_Key = pd.merge(synbio_one_col, chassis_one_col, on='InChI-Key', how='inner').reset_index(drop=True)
    print('I have merged')
    # Returns the InChI Key names by referencing the index
    # Converts array into a single column list
    # Removes any spaces
    shared_InChI_Key['InChI-Key'] = shared_InChI_Key['InChI-Key'].str.strip()
    # Finds unique InChI Keys in the list
    unique_chassis_InChI_Key = shared_InChI_Key['InChI-Key'].drop_duplicates()
    unique_chassis_InChI_Key.reset_index(drop=True)
    # Saves list of InChI Keys
    unique_chassis_InChI_Key = pd.DataFrame(unique_chassis_InChI_Key, columns=['InChI-Key'])
    # Removes the common InChI-Keys such as proton, ATP, and saves the list
    chassis_inchi_keys_translated = relevant_compounds(unique_chassis_InChI_Key)
    chassis_inchi_keys_translated.to_csv(sb_name + 'mutualism1.txt', header=True, index=True, sep='\t')
    return chassis_inchi_keys_translated


def mutualism2_modified_pathway(sb_name, pathway, chassis_bsm, synbio_bsm):
    # Opens MetaCyc list of all reactions
    metacyc_all_rxns = pd.read_csv(pathway + '/All-reactions-of-MetaCyc.txt',
                                   delimiter='\t', header=0, index_col=0)
    # Converts into database
    metacyc_all_rxns = pd.DataFrame(metacyc_all_rxns)
    print('Old column names', metacyc_all_rxns.columns)
    metacyc_all_rxns.columns = ['EC-Number', 'Substrates', 'Substrates InChI-Key', 'Reactants',
                                'Reactants InChI-Key', 'Products', 'Products InChI-Key']
    # String processing of dataframe by removing 'EC-' from the start of the EC number
    print('Here is the reactants column', metacyc_all_rxns['Reactants InChI-Key'])
    metacyc_all_rxns['EC-Number'] = metacyc_all_rxns['EC-Number'].str.replace('EC-', '', regex=False)
    # Merges based on EC number to create a list of reactions/substrates occurring in top match
    top_match_rxns = pd.merge(chassis_bsm, metacyc_all_rxns, on='EC-Number', how='inner')
    print(top_match_rxns.columns)
    print(top_match_rxns['Reactants InChI-Key'])
    # Turn on to save the list of substrates found in top match
    top_match_rxns.to_csv('topmatch_' + sb_name + '_chassis_all_rxns.txt', header=True, index=True, sep='\t')
    # Isolates the InChI Key column and splits the column based on // to isolate all substrates for top match
    chassis_InChI_Key = pd.DataFrame(top_match_rxns['Reactants InChI-Key'].astype(str).str.split('//', expand=True))
    chassis_one_col = to_one_column(chassis_InChI_Key)
    chassis_one_col.to_csv(sb_name + '_chassis_all_rxns.txt', header=True, index=True, sep='\t')
    # Merges based on EC number to create a list of reactions/substrates occurring in synbio
    synbio_rxns = pd.merge(synbio_bsm, metacyc_all_rxns, on='EC-Number', how='inner')
    # Merges the synbio binary summary matrix with the metacyc list to find the InChI-Key lists
    # Splits InChI-Keys based on the '//' separator
    synbio_InChI_Key = pd.DataFrame(synbio_rxns['Products InChI-Key'].astype(str).str.split('//', expand=True))
    synbio_one_col = to_one_column(synbio_InChI_Key)
    # Turn on to save the list of substrates found in synbio
    # synbio_one_col.to_csv(sb_name+'_synbio_all_rxns.txt', header=True, index= True, sep='\t')
    # Isolates the InChI Key column and splits the column based on // to isolate all substrates for synbio
    # Creates an array of InChI Keys of substrates that can be found in both organisms, saves the array as index
    shared_InChI_Key = pd.merge(synbio_one_col, chassis_one_col, on='InChI-Key', how='inner').reset_index(drop=True)
    print('I have merged')
    # Returns the InChI Key names by referencing the index
    # Converts array into a single column list
    # Removes any spaces
    shared_InChI_Key['InChI-Key'] = shared_InChI_Key['InChI-Key'].str.strip()
    # Finds unique InChI Keys in the list
    unique_chassis_InChI_Key = shared_InChI_Key['InChI-Key'].drop_duplicates()
    unique_chassis_InChI_Key.reset_index(drop=True)
    # Saves list of InChI Keys
    unique_chassis_InChI_Key = pd.DataFrame(unique_chassis_InChI_Key, columns=['InChI-Key'])
    # Removes the common InChI-Keys such as proton, ATP, and saves the list
    chassis_inchi_keys_translated = relevant_compounds(unique_chassis_InChI_Key)
    chassis_inchi_keys_translated.to_csv(sb_name + '_mutualism2.txt', header=True, index=True, sep='\t')
    # unique_chassis_InChI_Key.to_csv(sb_name+'_synbiovschassis_inchikey.txt', header=True, index= True, sep='\t')
    print('Top Match vs. Synbio InChI Key Substrates Analysis Is Complete')
    return chassis_inchi_keys_translated


# Inputs for this function include: dataframe which lists all the different EC numbers found in ec_locator(),
# the location where the CompetitorFind framework is saved
# Output for this function is dataframe which contain InChIKeys of substrates that are shared in the modified pathways
# of the synbio and top match, output for one is input for another = mutualism
# This method is not effective for unilateral changes (e.i. synbio is created by turning on EC numbers from top match)
# Method is preferred if synbio and top match both have EC numbers turned on and off

def competition_modified_pathway(sb_name, pathway, different_ECs):
    # Opens MetaCyc list of all reactions
    metacyc_all_rxns = pd.read_csv(pathway + '/All-reactions-of-MetaCyc.txt',
                                   delimiter='\t', header=0, index_col=0)
    # Converts into database
    metacyc_all_rxns = pd.DataFrame(metacyc_all_rxns)
    metacyc_all_rxns.columns = ['EC-Number', 'Substrates', 'Substrates InChI-Key', 'Reactants', 'Reactants InChI-Key',
                                'Products', 'Products InChI-Key']
    metacyc_all_rxns['EC-Number'] = metacyc_all_rxns['EC-Number'].str.replace('EC-', '', regex=False)
    # Finds the reactions associated with the modified pathway
    merged_rxns = pd.merge(different_ECs, metacyc_all_rxns, on='EC-Number', how='inner')
    # Turn on to save the list of all reactions that occur due to modified pathways
    # merged_rxns.to_csv(sb_name+'_altered_pathway_merged_rxns.txt', header=True, index=True, sep='\t')
    # Splits dataframes based on whether the top match EC numbers are being turned on
    top_match_modification = merged_rxns[merged_rxns.iloc[:, 1] == 1]
    # Isolates the InChI Key column and splits the column based on // to isolate all substrates for top match
    top_match_mod_InChIKey = top_match_modification['Reactants InChI-Key'].astype(str).str.split('//', expand=True)
    top_match_one_col = to_one_column(top_match_mod_InChIKey)
    # Splits dataframes based on whether the synbio EC Numbers are being turned on
    synbio_modification = merged_rxns[merged_rxns.iloc[:, 2] == 1]
    # Isolates the InChI Key column and splits the column based on // to isolate all substrates for synbio
    synbio_mod_InChIKey = synbio_modification['Reactants InChI-Key'].astype(str).str.split('//', expand=True)
    synbio_mod_one_col = to_one_column(synbio_mod_InChIKey)
    # Finds shared InChI Keys in the modified pathways, saves array as index of occurrence
    # Returns the InChI Key names by referencing the index
    pathway_shared_InChI_Key = pd.merge(top_match_one_col, synbio_mod_one_col, on=['InChI-Key'],
                                        how='inner').reset_index(drop=True)
    pathway_shared_InChI_Key.to_csv('BEFORE_EDITTING.txt', header=True, index=True, sep='\t')
    # Removes white space from the column
    pathway_shared_InChI_Key['InChI-Key'] = pathway_shared_InChI_Key['InChI-Key'].str.strip()
    # Finds unique InChI Keys in the list
    unique_path_InChI_Key = pathway_shared_InChI_Key['InChI-Key'].drop_duplicates()
    unique_path_InChI_Key.reset_index(drop=True)
    unique_path_InChI_Key = pd.DataFrame(unique_path_InChI_Key, columns=['InChI-Key'])
    # Saves list of InChI Keys
    unique_pathway_translated = relevant_compounds(unique_path_InChI_Key)
    unique_pathway_translated.to_csv(sb_name + '_modifiedpathway_inchikey.txt', header=True, index=True, sep='\t')
    print('Modified Pathway Substrates Analysis Is Complete')
    return unique_path_InChI_Key


##====================================================================================================================##

pathway = '/home/anna/Desktop/EcoGenoRisk/HazID/CompetitorFind'
# synbio = 'Aquificota_Actinobacteria_Chimera'
synbio = 'Chimera1'
path = '/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/TestCases/Chimera1'
# path = '/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/Similar'
# synbio = 'E_Coli_Chimera'
# synbio = 'E_Coli_Chimera_1172_Off'
# Changes directory to local
os.chdir(pathway)
# Sends to ec_locator() which finds the EC BSM for top match and synbio. Returns the BSM for both organisms
[different_ECs, chassis_bsm, synbio_bsm, top_match] = ec_locator(synbio, path)
# Sends to function to compare shared InChiKeys in synbio and top match
# Resulting matrix should be 1000-10,000 InChiKeys as this method is the "catch-all"
# Returns the shared InChiKeys without white space or duplicates
individual_genome_rxns = substrate_changes_synbio_v_chassis(synbio, pathway, synbio_bsm, chassis_bsm)
# Converts InChiKey notation to laymen terms
translated_individual_rxns = inchikey_to_conventional_names(individual_genome_rxns)
# Saves output
translated_individual_rxns.to_csv(synbio + '_names.txt', header=True, index=True, sep='\t')
# Sends EC changes between synbio and top match to find shared InChiKeys
# Resulting matrix should be 10-100 InChiKeys as method is more precise looks at InChiKeys within the modified pathway
# Returns the shared InChiKeys without white space or duplicates
mutualism1 = mutualism1_modified_pathway(synbio, pathway, chassis_bsm,synbio_bsm)
mutualism2 = mutualism2_modified_pathway(synbio, pathway, chassis_bsm,synbio_bsm)
pathway_modification = competition_modified_pathway(synbio, pathway, different_ECs)
# Converts InChiKey notation to laymen terms
translated_pathway = inchikey_to_conventional_names(pathway_modification)
# Saves output
translated_pathway.to_csv(synbio + '_modifiedpathway_names.txt', header=True, index=True, sep='\t')
# Finds InChiKeys that are present in both sets
non_repeated_InChIKeys = set(individual_genome_rxns) ^ set(pathway_modification)
# Isolates only the column with the listed InChI-Keys
non_repeated_InChIKeys = pd.DataFrame(non_repeated_InChIKeys, columns=['InChI-Key'])
# Translates from InChI-Key format to conventional names
duplicitous_compounds = inchikey_to_conventional_names(non_repeated_InChIKeys)
# Saves as a CSV file
duplicitous_compounds.to_csv(synbio + '_duplicitous_compounds_translated.txt', header=True, index=True, sep='\t')
# Sends to scoring() function to render score. Mutualism is substracted from the overall score as it promotes microbial
# growth and sustains biodiversity. Weighting was decided based on number of compounds produced per run
score = scoring(individual_genome_rxns, pathway_modification, non_repeated_InChIKeys, mutualism1, mutualism2)
# distance = finding_distance(top_match)
data = {'Individual Reactions': individual_genome_rxns.shape[0], 'Modified Pathway': pathway_modification.shape[0],
        'Set Overlap': non_repeated_InChIKeys.shape[0], 'Mutualism 1 Observed': mutualism1.shape[0],
        'Mutualism 2 Observed': mutualism2.shape[0], 'Overall Score': score}  # , 'Distance': distance}

overview = pd.DataFrame(data, index=[0])
overview.to_csv(synbio+'_CompetitorFind_Analysis_invasive.txt', header=True, index=True, sep='\t')
