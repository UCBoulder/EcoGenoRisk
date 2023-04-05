import numpy
import pandas as pd
import numpy as np
import os as os
from scoring import scoring


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
    # Removes row if InChIKEy for water is present
    df.drop(index_water, inplace=True)
    print(df)
    print('water removed')
    # Finds the indexes where InChiKey for NADPH are present
    index_NADPH = df[df['InChI-Key'] == 'InChIKey=ACFIXJIJDZMPPO-NNYOXOHSSA-J'].index
    # Removes row if InChIKey for NADPH is present
    df.drop(index_NADPH, inplace=True)
    print(df)
    print('NADPH removed')
    # Finds the indexes where InChiKey for NADP are present
    index_NADP = df[df['InChI-Key'] == 'InChIKey=XJLXINKUBYWONI-NNYOXOHSSA-K'].index
    # Removes row if InChIKey for NADP is present
    df.drop(index_NADP, inplace=True)
    print(df)
    print('NADP removed')
    # Finds the indexes where InChiKey for Proton are present
    index_proton = df[df['InChI-Key'] == 'InChIKey=GPRLSGONYQIRFK-UHFFFAOYSA-N'].index
    # Removes row if InChIKey for Proton is present
    df.drop(index_proton, inplace=True)
    print(df)
    print('proton removed')
    # Resets index
    #df.reset_index(drop=True)
    # Returns list of relevant InChiKeys
    df2 = df
    return df2
##====================================================================================================================##

# Translates InChI-Keys into conventional compound names
# Returns dataframe of Inchikeys and its translation
def inchikey_to_conventional_names(df):
    key = pd.read_csv('InchiKeystoCompoundNames.txt', delimiter='\t', header=0, index_col=None)
    translated_inchikeys = pd.merge(df, key, on='InChI-Key', how='inner')
    return translated_inchikeys

##====================================================================================================================##

# Inputs for this function include: location where the CompetitorFind framework is saved,
# top matches original EC BSM, synbio EC BSM
# Output for this function is dataframe which contain InChIKeys of substrates that are shared between the synbio
# organism and the top match
# Method is a "catch all", no matter what kind of changes occur in the EC BSM
def substrate_changes_synbio_v_chassis(sb_name, pathway, chassis_bsm, synbio_bsm):
    # Opens MetaCyc list of all reactions
    print('Starting overall substrate overlap analysis')
    metacyc_all_rxns = pd.read_csv(pathway + '/All-reactions-of-MetaCyc.txt',
                                   delimiter='\t', header=0, index_col=0)
    # Converts into database
    metacyc_all_rxns = pd.DataFrame(metacyc_all_rxns)
    metacyc_all_rxns.columns=['EC-Number', 'Substrates','Substrates InChI-Key', 'Reactants',
                                                               'Reactants InChI-Key','Products', 'Products InChI-Key']
    # String processing of dataframe by removing 'EC-' from the start of the EC number
    print('Here is the reactants column',metacyc_all_rxns['Reactants InChI-Key'])
    metacyc_all_rxns['EC-Number'] = metacyc_all_rxns['EC-Number'].str.replace('EC-', '', regex=False)
    # Merges based on EC number to create a list of reactions/substrates occurring in top match
    top_match_rxns = pd.merge(chassis_bsm, metacyc_all_rxns, on='EC-Number', how='inner')
    print(top_match_rxns.columns)
    print(top_match_rxns['Reactants InChI-Key'])
    # Turn on to save the list of substrates found in top match
    #top_match_rxns.to_csv('topmatch_'+sb_name+'_chassis_all_rxns.txt', header=True, index=True, sep='\t')
    # Isolates the InChI Key column and splits the column based on // to isolate all substrates for top match
    chassis_InChI_Key = pd.DataFrame(top_match_rxns['Reactants InChI-Key'].astype(str).str.split(' // ', expand=True))
    chassis_one_col = to_one_column(chassis_InChI_Key)
    #chassis_one_col.to_csv(sb_name+'_chassis_all_rxns.txt', header=True, index= True, sep='\t')
    # Merges based on EC number to create a list of reactions/substrates occurring in synbio

    synbio_rxns = pd.merge(synbio_bsm, metacyc_all_rxns, on='EC-Number', how='inner')
    synbio_InChI_Key = pd.DataFrame(synbio_rxns['Reactants InChI-Key'].astype(str).str.split(' // ', expand=True))
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
    print('Overall Reactions', shared_InChI_Key.shape)
    #shared_InChI_Key.to_csv('overall_reactions.txt', sep='\t')
    # For optimization, remove the blank rows first before further filtering and refining compound list
    #shared_InChI_Key['InChI-Key']=shared_InChI_Key[[shared_InChI_Key['InChI-Key'].str.contains('')==False]]
    print('Overall Reactions', shared_InChI_Key.shape)
    shared_InChI_Key['InChI-Key'] = shared_InChI_Key['InChI-Key'].str.strip()
    #shared_InChI_Key.to_csv('stripped_white_space.txt',sep='\t')
    # Finds unique InChI Keys in the list
    unique_chassis_InChI_Key = shared_InChI_Key['InChI-Key'].drop_duplicates()
    #unique_chassis_InChI_Key.to_csv('dropped_duplicated.txt',sep='\t')
    unique_chassis_InChI_Key.reset_index(drop=True)
    # Saves list of InChI Keys
    unique_chassis_InChI_Key = pd.DataFrame(unique_chassis_InChI_Key, columns=['InChI-Key'])
    unique_chassis_InChI_Key = unique_chassis_InChI_Key[unique_chassis_InChI_Key['InChI-Key'].str.contains('None')==False]
    #unique_chassis_InChI_Key.to_csv('removed_None.txt',sep='\t')
    unique_chassis_InChI_Key = unique_chassis_InChI_Key[unique_chassis_InChI_Key['InChI-Key'].str.contains('na')==False]
    #unique_chassis_InChI_Key.to_csv('removed_na.txt', sep='\t')
    chassis_inchi_keys_translated = relevant_compounds(unique_chassis_InChI_Key)
    print(chassis_inchi_keys_translated)
    #chassis_inchi_keys_translated.to_csv(sb_name + '_synbiovschassis_inchikey.txt', header=True, index=True, sep='\t')
    chassis_inchi_keys_translated.to_csv(sb_name+'_synbiovschassis_inchikey_new.txt', header=True, index= True, sep='\t')
    print('Top Match vs. Synbio InChI Key Substrates Analysis Is Complete')
    return unique_chassis_InChI_Key


##====================================================================================================================##

# Inputs for this function include: dataframe which lists all the different EC numbers found in ec_locator(),
# the location where the CompetitorFind framework is saved
# Output for this function is dataframe which contain InChIKeys of substrates that are shared in the modified pathways
# of the synbio and top match
# This method is not effective for unilateral changes (e.i. synbio is created by turning on EC numbers from top match)
# Method is preferred if synbio and top match both have EC numbers turned on and off

def competition_modified_pathway(sb_name, pathway, different_ECs):
    # Opens MetaCyc list of all reactions
    metacyc_all_rxns = pd.read_csv(pathway + '/All-reactions-of-MetaCyc.txt',
                                   delimiter='\t', header=0, index_col=0)
    # Converts into database
    metacyc_all_rxns = pd.DataFrame(metacyc_all_rxns)
    metacyc_all_rxns.columns=['EC-Number', 'Substrates','Substrates InChI-Key','Reactants','Reactants InChI-Key',
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
    # Removes white space from the column
    pathway_shared_InChI_Key['InChI-Key'] = pathway_shared_InChI_Key['InChI-Key'].str.strip()
    # Finds unique InChI Keys in the list
    unique_path_InChI_Key = pathway_shared_InChI_Key['InChI-Key'].drop_duplicates()
    unique_path_InChI_Key.reset_index(drop=True)
    unique_path_InChI_Key = pd.DataFrame(unique_path_InChI_Key, columns=['InChI-Key'])
    unique_path_InChI_Key = unique_path_InChI_Key.dropna()
    unique_path_InChI_Key =unique_path_InChI_Key[unique_path_InChI_Key['InChI-Key'].str.contains('None')==False]
    unique_path_InChI_Key =unique_path_InChI_Key[unique_path_InChI_Key['InChI-Key'].str.contains('')==False]
    unique_path_InChI_Key =unique_path_InChI_Key[unique_path_InChI_Key['InChI-Key'].str.contains(' ')==False]
    unique_path_InChI_Key =unique_path_InChI_Key[unique_path_InChI_Key['InChI-Key'].str.contains('na')==False]
    unique_path_InChI_Key =unique_path_InChI_Key[unique_path_InChI_Key['InChI-Key'].str.contains('NaN')==False]
    unique_path_InChI_Key =unique_path_InChI_Key[unique_path_InChI_Key['InChI-Key'].str.contains('NAN')==False]
    # Saves list of InChI Keys
    #unique_pathway_translated = relevant_compounds(unique_path_InChI_Key)
    #unique_pathway_translated.to_csv(sb_name + '_modifiedpathway_inchikey.txt', header=True, index=True, sep='\t')
    print('Modified Pathway Substrates Analysis Is Complete')
    return unique_path_InChI_Key

##====================================================================================================================##

def mutualism1_modified_pathway(sb_name, pathway, chassis_bsm, synbio_bsm):
    # Opens MetaCyc list of all reactions
    metacyc_all_rxns = pd.read_csv(pathway + '/All-reactions-of-MetaCyc.txt',
                                   delimiter='\t', header=0, index_col=0)
    # Converts into database
    metacyc_all_rxns = pd.DataFrame(metacyc_all_rxns)
    print('Old column names', metacyc_all_rxns.columns)
    metacyc_all_rxns.columns=['EC-Number', 'Substrates','Substrates InChI-Key', 'Reactants',
                                                               'Reactants InChI-Key','Products', 'Products InChI-Key']
    # String processing of dataframe by removing 'EC-' from the start of the EC number
    print('Here is the reactants column',metacyc_all_rxns['Reactants InChI-Key'])
    metacyc_all_rxns['EC-Number'] = metacyc_all_rxns['EC-Number'].str.replace('EC-', '', regex=False)
    # Merges based on EC number to create a list of reactions/substrates occurring in top match
    top_match_rxns = pd.merge(chassis_bsm, metacyc_all_rxns, on='EC-Number', how='inner')
    print(top_match_rxns.columns)
    print(top_match_rxns['Reactants InChI-Key'])
    # Turn on to save the list of substrates found in top match
    #top_match_rxns.to_csv('topmatch_'+sb_name+'_chassis_all_rxns.txt', header=True, index=True, sep='\t')
    # Isolates the InChI Key column and splits the column based on // to isolate all substrates for top match
    chassis_InChI_Key = pd.DataFrame(top_match_rxns['Reactants InChI-Key'].astype(str).str.split('//', expand=True))
    chassis_one_col = to_one_column(chassis_InChI_Key)
    #chassis_one_col.to_csv(sb_name+'_chassis_all_rxns.txt', header=True, index= True, sep='\t')
    # Merges based on EC number to create a list of reactions/substrates occurring in synbio

    synbio_rxns = pd.merge(synbio_bsm, metacyc_all_rxns, on='EC-Number', how='inner')
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
    unique_chassis_InChI_Key = unique_chassis_InChI_Key[unique_chassis_InChI_Key['InChI-Key'].str.contains('None')==False]
    unique_chassis_InChI_Key = unique_chassis_InChI_Key[unique_chassis_InChI_Key['InChI-Key'].str.contains('na')==False]

    #chassis_inchi_keys_translated = relevant_compounds(unique_chassis_InChI_Key)
    #chassis_inchi_keys_translated.to_csv(sb_name + '_synbiovschassis_inchikey.txt', header=True, index=True, sep='\t')
    unique_chassis_InChI_Key.to_csv(sb_name+'mutualism1.txt', header=True, index= True, sep='\t')
    print('Top Match vs. Synbio InChI Key Substrates Analysis Is Complete')
    return unique_chassis_InChI_Key

##====================================================================================================================##

def mutualism2_modified_pathway(sb_name, pathway, chassis_bsm, synbio_bsm):
    # Opens MetaCyc list of all reactions
    metacyc_all_rxns = pd.read_csv(pathway + '/All-reactions-of-MetaCyc.txt',
                                   delimiter='\t', header=0, index_col=0)
    # Converts into database
    metacyc_all_rxns = pd.DataFrame(metacyc_all_rxns)
    print('Old column names', metacyc_all_rxns.columns)
    metacyc_all_rxns.columns=['EC-Number', 'Substrates','Substrates InChI-Key', 'Reactants',
                                                               'Reactants InChI-Key','Products', 'Products InChI-Key']
    # String processing of dataframe by removing 'EC-' from the start of the EC number
    print('Here is the reactants column',metacyc_all_rxns['Reactants InChI-Key'])
    metacyc_all_rxns['EC-Number'] = metacyc_all_rxns['EC-Number'].str.replace('EC-', '', regex=False)
    # Merges based on EC number to create a list of reactions/substrates occurring in top match
    top_match_rxns = pd.merge(chassis_bsm, metacyc_all_rxns, on='EC-Number', how='inner')
    print(top_match_rxns.columns)
    print(top_match_rxns['Reactants InChI-Key'])
    # Turn on to save the list of substrates found in top match
    #top_match_rxns.to_csv('topmatch_'+sb_name+'_chassis_all_rxns.txt', header=True, index=True, sep='\t')
    # Isolates the InChI Key column and splits the column based on // to isolate all substrates for top match
    chassis_InChI_Key = pd.DataFrame(top_match_rxns['Reactants InChI-Key'].astype(str).str.split('//', expand=True))
    chassis_one_col = to_one_column(chassis_InChI_Key)
    #chassis_one_col.to_csv(sb_name+'_chassis_all_rxns.txt', header=True, index= True, sep='\t')
    # Merges based on EC number to create a list of reactions/substrates occurring in synbio

    synbio_rxns = pd.merge(synbio_bsm, metacyc_all_rxns, on='EC-Number', how='inner')
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
    unique_chassis_InChI_Key = unique_chassis_InChI_Key[unique_chassis_InChI_Key['InChI-Key'].str.contains('None')==False]
    unique_chassis_InChI_Key = unique_chassis_InChI_Key[unique_chassis_InChI_Key['InChI-Key'].str.contains('na')==False]

    #chassis_inchi_keys_translated = relevant_compounds(unique_chassis_InChI_Key)
    #chassis_inchi_keys_translated.to_csv(sb_name + '_synbiovschassis_inchikey.txt', header=True, index=True, sep='\t')
    unique_chassis_InChI_Key.to_csv(sb_name+'mutualism2.txt', header=True, index= True, sep='\t')
    print('Top Match vs. Synbio InChI Key Substrates Analysis Is Complete')
    return unique_chassis_InChI_Key