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
    print(top_match_name, '\n', top_match_bsm)
    print(sb_name, '\n', synbio_bsm)
    outer_merge = pd.merge(top_match_bsm,synbio_bsm, on='EC-Number', how ='outer')
    outer_merge = outer_merge[outer_merge[sb_name]!=outer_merge[top_match_name]]
    #outer_merge.to_csv('outer_merge', header=True, index=True, sep='\t')
    return outer_merge

def changes_in_rxn(pathway, outer_merge):
    metacyc_all_rxns = pd.read_csv(pathway+'/20220802_Copy-of-All-reactions-inhibitors-of-MetaCyc.txt',
                                                delimiter='\t', header = 0, index_col = 0)
    metacyc_all_rxns = pd.DataFrame(metacyc_all_rxns)
    print(metacyc_all_rxns)
    print(type(metacyc_all_rxns))
    metacyc_all_rxns['EC-Number'] = metacyc_all_rxns['EC-Number'].str.replace('EC-', '', regex = False)
    merged_rxns = pd.merge(outer_merge, metacyc_all_rxns, on = 'EC-Number', how = 'inner')
    merged_rxns.to_csv('merged_rxns', header=True, index= True, sep='\t')
    return merged_rxns['InChI-Key Inhibitor	InChI Inhibitor']

outer_merge = ec_locator('E_Coli_Chimeria', '/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/Similar')
pathway = '/home/anna/Desktop/EcoGenoRisk/HazID/CompetitorFind'
changes_in_rxn(pathway, outer_merge)