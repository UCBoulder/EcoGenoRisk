#_______________________________________________________________________________
# Checks for matches between two functional profiles, or the enzymatic overlap between two organisms

import pandas as pd



kleb = pd.read_csv('/home/anna/Downloads/kleb_funct_profile', delimiter = ' ')
klbe2 = kleb.T
# print(klbe2)
new_header = klbe2.iloc[0]
kleb3 = klbe2[1:]
kleb3.columns=new_header
print('kleb3 \n', kleb3)


Paen = pd.read_csv('/home/anna/Downloads/Paenibacllius_functional_profile.txt', delimiter = ' ')
paen2 = Paen.T
new_header2 = paen2.iloc[0]
paen3 = paen2[1:]
paen3.columns=new_header2
pane4 = paen3.rename(columns={'3300031908_4_matches.tsv': 'panticola'})
print(paen3)

combined = pd.concat([kleb3, pane4], axis = 1)
print(combined)

# combined.to_excel('Combined.xlsx')
# reset = combined.reset_index
# print(reset)
print(len(combined.query('`Synbio_matches.tsv` == panticola')))

present = combined[(combined['Synbio_matches.tsv'] !=0) & (combined['panticola'] !=0)]
print(present)