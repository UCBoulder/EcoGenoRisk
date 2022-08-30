import numpy
import pandas as pd
import numpy as np

def ec_locator(sb_name, pathway):
    find_top_match = pd.read_csv(pathway+'/Difference_Based_Comparison_Score.txt', delimiter='\t', header=0, index_col=0)
    top_match_row = find_top_match.head(1).reset_index(drop=True)
    top_match_name = top_match_row.iloc[0,0]
    ec_bsm = pd.read_csv(pathway+'/complete_binary_matrix.txt', delimiter='\t', header=0, index_col=0)
    top_match_bsm = pd.DataFrame(ec_bsm.loc[top_match_name], columns= [top_match_name])
    top_match_bsm.index.name =  'EC-Number'
    synbio_bsm = pd.DataFrame(ec_bsm.loc[sb_name], columns=[sb_name])
    synbio_bsm.index.name = 'EC-Number'
    outer_merge = pd.merge(top_match_bsm,synbio_bsm, on='EC-Number', how ='outer')
    outer_merge = outer_merge[outer_merge[sb_name]!=outer_merge[top_match_name]]
    #outer_merge.to_csv('outer_merge', header=True, index=True, sep='\t')
    return outer_merge, top_match_bsm, top_match_name, sb_name

def changes_in_rxn(pathway, outer_merge, chassis_bsm):
    metacyc_all_rxns = pd.read_csv(pathway+'/All-reactions-of-MetaCyc.txt',
                                                delimiter='\t', header = 0, index_col = 0)
    metacyc_all_rxns = pd.DataFrame(metacyc_all_rxns)
    metacyc_all_rxns['EC-Number'] = metacyc_all_rxns['EC-Number'].str.replace('EC-', '', regex = False)
    top_match_rxns = pd.merge(chassis_bsm, metacyc_all_rxns, on = 'EC-Number', how = 'inner')
    #top_match_rxns.to_csv('chassis_all_rxns.txt', header=True, index= True, sep='\t')
    merged_rxns = pd.merge(outer_merge, metacyc_all_rxns, on = 'EC-Number', how = 'inner')
    #merged_rxns.to_csv('altered_pathway_merged_rxns.txt', header=True, index= True, sep='\t')

    # Compares chassis to modified pathway
    chassis_InChI_Key = pd.DataFrame(top_match_rxns['InChI-Key'].str.split('//', expand = True))
    # chassis_InChI_Key.to_csv('parent_inchikey.txt', header=False, index = False, sep = '\t')
    altered_pathway_InChI_Key = pd.DataFrame(merged_rxns['InChI-Key'].str.split('//', expand = True))
    # altered_pathway_InChI_Key.to_csv('modified_inchikey.txt', header=False, index = False, sep = '\t')
    chassis_vs_mod = np.isin(altered_pathway_InChI_Key, chassis_InChI_Key)
    chassis_shared_InChI_Key = altered_pathway_InChI_Key[chassis_vs_mod]
    ###
    chassis_vs_synbio = []
    for column in chassis_shared_InChI_Key:
        chassis_synbio_single_col = chassis_shared_InChI_Key[column].tolist()
        chassis_vs_synbio += chassis_synbio_single_col
    chassis_one_column_df = pd.DataFrame(chassis_vs_synbio, columns = ['InChI_Key'])
    unique_chassis_InChI_Key = chassis_one_column_df['InChI_Key'].unique()
    unique_chassis_InChI_Key = unique_chassis_InChI_Key[unique_chassis_InChI_Key!='InChIKey=XLYOFNOQVPJJNP-UHFFFAOYSA-N']
    unique_chassis_InChI_Key = unique_chassis_InChI_Key[unique_chassis_InChI_Key!='None']
    numpy.savetxt('synbiovschassis_inchikey.txt', unique_chassis_InChI_Key, fmt='%s', delimiter = '\t')
    # Compares modified pathway
    top_match_modification  = merged_rxns[merged_rxns.iloc[:,1].values==1]
    print(top_match_modification)
    top_match_mod_InChIKey =  top_match_modification['InChI-Key'].str.split('//', expand = True)
    synbio_modification = merged_rxns[merged_rxns.iloc[:,2].values==1]
    synbio_mod_InChIKey = synbio_modification['InChI-Key'].str.split('//', expand = True)
    mod_pathway = np.isin(top_match_mod_InChIKey,synbio_mod_InChIKey)
    pathway_shared_InChI_Key = top_match_mod_InChIKey[mod_pathway]

    pathway_comp = []
    for column in pathway_shared_InChI_Key:
        pathway_list = pathway_shared_InChI_Key[column].tolist()
        pathway_comp += pathway_list
    modified_pathway_column = pd.DataFrame(pathway_comp, columns=['InChI_Key'])
    unique_path_InChI_Key = modified_pathway_column['InChI_Key'].unique()
    unique_path_InChI_Key = unique_path_InChI_Key[
        unique_path_InChI_Key != 'InChIKey=XLYOFNOQVPJJNP-UHFFFAOYSA-N']
    unique_path_InChI_Key = unique_path_InChI_Key[unique_path_InChI_Key != 'None']
    numpy.savetxt('modifiedpathway_inchikey.txt', unique_path_InChI_Key, fmt='%s', delimiter='\t')
    return unique_chassis_InChI_Key, unique_path_InChI_Key


[outer_merge, chassis_bsm, chassis, synbio] = ec_locator('E_Coli_Chimeria', '/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/Similar')
pathway = '/home/anna/Desktop/EcoGenoRisk/HazID/CompetitorFind'
[individual_genome_rxns, pathway_modifcation]= changes_in_rxn(pathway, outer_merge, chassis_bsm)
print(individual_genome_rxns, '\n', pathway_modifcation)