import numpy
import pandas as pd
import numpy as np


##====================================================================================================================##

# Inputs for this function are the pathway where synbio analysis documents were saved, and the name of the synbio organism
# Outputs include: a data frame which contains the EC numbers that are variable between the two genomes and which
# organism the EC numbers occur, the top match's genome name and EC BSM
# This function tracks changes in the EC BSM between the synbio organism and the top match found in the function
# genome_to_genome_diffcomp() from the _5_synbio_distance script
def ec_locator(sb_name, pathway):
    # Reads the document that scores all genomes based on EC expression
    find_top_match = pd.read_csv(pathway + '/Difference_Based_Comparison_Score.txt', delimiter='\t', header=0,
                                 index_col=0)
    # Finds the top match row from the document. This organism is functionally most similar to the synbio organism
    top_match_row = find_top_match.head(1).reset_index(drop=True)
    # Extracts the top match genome name
    top_match_name = top_match_row.iloc[0, 0]
    print('Top Match Found is: ', top_match_name)
    # Reads document which has EC BSM for Bacteria, Archaea, and the synbio organism
    ec_bsm = pd.read_csv(pathway+'/complete_binary_matrix _CompFind.txt', delimiter='\t', header=0, index_col=0)
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
    # Finds occurances where the synbio and top match EC BSMs don't have the same binary values and extracts EC numbers
    different_ECs = different_ECs[different_ECs[sb_name] != different_ECs[top_match_name]]
    # Turn on to save the list of EC numbers that are different between the two genomes and where the EC number is present
    different_ECs.to_csv('different_EC_profile.txt', header=True, index=True, sep='\t')
    print('Locating Variable EC numbers Complete')
    return different_ECs, top_match_bsm, synbio_bsm


##====================================================================================================================##

# Inputs for this function include: location where the CompetitorFind framework is saved,
# top matches original EC BSM, synbio EC BSM
# Output for this function is dataframe which contain InChIKeys of substrates that are shared between the synbio
# organism and the top match
# Method is a "catch all", no matter what kind of changes occur in the EC BSM
def substrate_changes_synbio_v_chassis(pathway, chassis_bsm, synbio_bsm, different_ECs):
    # Opens MetaCyc list of all reactions
    metacyc_all_rxns = pd.read_csv(pathway + '/All-reactions-of-MetaCyc.txt',
                                   delimiter='\t', header=0, index_col=0)
    # Converts into database
    metacyc_all_rxns = pd.DataFrame(metacyc_all_rxns)
    # String processing of dataframe by removing 'EC-' from the start of the EC number
    metacyc_all_rxns['EC-Number'] = metacyc_all_rxns['EC-Number'].str.replace('EC-', '', regex=False)
    # Merges based on EC number to create a list of reactions/substrates occurring in top match
    top_match_rxns = pd.merge(chassis_bsm, metacyc_all_rxns, on='EC-Number', how='inner')
    # Turn on to save the list of substrates found in top match
    #top_match_rxns.to_csv('chassis_all_rxns.txt', header=True, index= True, sep='\t')
    # Isolates the InChI Key column and splits the column based on // to isolate all substrates for top match
    chassis_InChI_Key = pd.DataFrame(top_match_rxns['InChI-Key'].str.split('//', expand=True))
    # Merges based on EC number to create a list of reactions/substrates occurring in synbio
    synbio_rxns = pd.merge(synbio_bsm, metacyc_all_rxns, on='EC-Number', how='inner')
    # Turn on to save the list of substrates found in synbio
    #synbio_rxns.to_csv('synbio_all_rxns.txt', header=True, index= True, sep='\t')
    # Isolates the InChI Key column and splits the column based on // to isolate all substrates for synbio
    synbio_InChI_Key = pd.DataFrame(synbio_rxns['InChI-Key'].str.split('//', expand=True))
    # Creates an array of InChI Keys of substrates that can be found in both organisms, saves the array as index
    chassis_vs_syn = np.isin(synbio_InChI_Key, chassis_InChI_Key)
    # Returns the InChI Key names by referencing the index
    shared_InChI_Key = synbio_InChI_Key[chassis_vs_syn]
    # Converts array into a single column list
    chassis_vs_synbio = []
    for column in shared_InChI_Key:
        chassis_synbio_single_col = shared_InChI_Key[column].tolist()
        chassis_vs_synbio += chassis_synbio_single_col
    chassis_one_column_df = pd.DataFrame(chassis_vs_synbio, columns=['InChI_Key'])
    # Removes any spaces
    chassis_one_column_df['InChI_Key'] = chassis_one_column_df['InChI_Key'].str.strip()
    # Finds unique InChI Keys in the list
    unique_chassis_InChI_Key = chassis_one_column_df.drop_duplicates()
    # Excludes InChIKey for water to give other substrates a higher priority
    unique_chassis_InChI_Key = unique_chassis_InChI_Key[
        unique_chassis_InChI_Key['InChI_Key'].str.contains('InChIKey=XLYOFNOQVPJJNP-UHFFFAOYSA-N') == False]
    # Removes any InChI Keys that were listed as None
    unique_chassis_InChI_Key = unique_chassis_InChI_Key[
        unique_chassis_InChI_Key['InChI_Key'].str.contains('None') == False]
    # Saves list of InChI Keys
    numpy.savetxt('synbiovschassis_inchikey.txt', unique_chassis_InChI_Key, fmt='%s', delimiter='\t')
    print('Top Match vs. Synbio InChI Key Substrates Analysis Is Complete')
    return unique_chassis_InChI_Key


##====================================================================================================================##

# Inputs for this function include: dataframe which lists all the different EC numbers found in ec_locator(),
# the location where the CompetitorFind framework is saved
# Output for this function is dataframe which contain InChIKeys of substrates that are shared in the modified pathways
# of the synbio and top match
# This method is not effective for unilateral changes (e.i. synbio is created by turning on EC numbers from top match)
# Method is preferred if synbio and top match both have EC numbers turned on and off

def substrate_changes_modified_pathway(pathway, different_ECs):
    # Opens MetaCyc list of all reactions
    metacyc_all_rxns = pd.read_csv(pathway + '/All-reactions-of-MetaCyc.txt',
                                   delimiter='\t', header=0, index_col=0)
    # Converts into database
    metacyc_all_rxns = pd.DataFrame(metacyc_all_rxns)
    metacyc_all_rxns['EC-Number'] = metacyc_all_rxns['EC-Number'].str.replace('EC-', '', regex=False)
    # Finds the reactions associated with the modified pathway
    merged_rxns = pd.merge(different_ECs, metacyc_all_rxns, on='EC-Number', how='inner')
    # Turn on to save the list of all reactions that occur due to modified pathways
    merged_rxns.to_csv('altered_pathway_merged_rxns.txt', header=True, index = True, sep='\t')
    # Splits dataframes based on whether the top match EC numbers are being turned on
    top_match_modification = merged_rxns[merged_rxns.iloc[:, 1] == 1]
    # Isolates the InChI Key column and splits the column based on // to isolate all substrates for top match
    top_match_mod_InChIKey = top_match_modification['InChI-Key'].str.split('//', expand=True)
    # Splits dataframes based on whether the synbio EC Numbers are being turned on
    synbio_modification = merged_rxns[merged_rxns.iloc[:, 2] == 1]
    ## DELETE ONCE DONE TESTING, TRYING TO SEE IF THIS SHIT PROPERLY PARTITIONS THE DFS
    ##
    print('Top Match Partition: ', '\n', top_match_modification)
    print('Synbio Partition: ', '\n', synbio_modification)
    ##
    ##
    # Isolates the InChI Key column and splits the column based on // to isolate all substrates for synbio
    synbio_mod_InChIKey = synbio_modification['InChI-Key'].str.split('//', expand=True)
    # Finds shared InChI Keys in the modified pathways, saves array as index of occurrence
    mod_pathway = np.isin(top_match_mod_InChIKey, synbio_mod_InChIKey)
    # Returns the InChI Key names by referencing the index
    pathway_shared_InChI_Key = top_match_mod_InChIKey[mod_pathway]
    np.savetxt('raw_unprocessed_shared_keys.txt', pathway_shared_InChI_Key, fmt = '%s', delimiter='\t')
    # Converts array into a single column list
    pathway_comp = []
    for column in pathway_shared_InChI_Key:
        pathway_list = pathway_shared_InChI_Key[column].tolist()
        pathway_comp += pathway_list
    modified_pathway_column = pd.DataFrame(pathway_comp, columns=['InChI_Key'])
    print(modified_pathway_column)
    print(type(modified_pathway_column))
    # Removes white space from the column
    modified_pathway_column['InChI_Key'] = modified_pathway_column['InChI_Key'].str.strip()
    # Finds unique InChI Keys in the list
    unique_path_InChI_Key = modified_pathway_column.drop_duplicates()
    # Excludes InChIKey for water to give other substrates a higher priority
    unique_path_InChI_Key = unique_path_InChI_Key[
        unique_path_InChI_Key['InChI_Key'].str.contains('InChIKey=XLYOFNOQVPJJNP-UHFFFAOYSA-N') == False]
    # Removes any InChI Keys that were listed as None
    unique_path_InChI_Key = unique_path_InChI_Key[
        unique_path_InChI_Key['InChI_Key'].str.contains('None') == False]
    # Saves list of InChI Keys
    numpy.savetxt('modifiedpathway_inchikey.txt', unique_path_InChI_Key, fmt='%s', delimiter='\t')
    print('Modified Pathway Substrates Analysis Is Complete')
    return unique_path_InChI_Key
##====================================================================================================================##

#synbio = 'Aquificota_Actinobacteria_Chimera'
#synbio = 'Aquificota_Actinobacteria_Chimera_935_Off'
#path = '/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/Different'
path = '/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/Similar'
#synbio = 'E_Coli_Chimera'
synbio = 'E_Coli_Chimera_1172_Off'
[different_ECs, chassis_bsm, synbio_bsm] = ec_locator(synbio, path)
pathway = '/home/anna/Desktop/EcoGenoRisk/HazID/CompetitorFind'
individual_genome_rxns = substrate_changes_synbio_v_chassis(pathway, synbio_bsm, chassis_bsm, different_ECs)
pathway_modifcation = substrate_changes_modified_pathway(pathway, different_ECs)
print(individual_genome_rxns, '\n', pathway_modifcation)
