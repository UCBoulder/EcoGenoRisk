
import pandas as pd
import numpy as np
import os as os

def scoring(df1, df2, df3, df4, df5):
    overall_sim = df1.size #gives the rows*columns
    print('size of overall_sim is \n', overall_sim)
    pathway_sim = df2.size
    common_sim = df3.size
    mutualism = df4.size
    mutualism2 = df5.size
    score = overall_sim/10000 + pathway_sim/100 + common_sim/10 - mutualism/200 - mutualism2/200
    return score

def finding_distance(genome2):
    complete_binary = pd.read_csv(
        '/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/dm_df_labelled.txt',
        delimiter='\t', header=0)
    complete_binary.columns = ['Genome', 'Distance']
    complete_binary['Genome'] = complete_binary['Genome'].str.strip('(')
    complete_binary['Genome'] = complete_binary['Genome'].str.strip(')')
    complete_binary['Genome'] = complete_binary['Genome'].str.strip(',')
    complete_binary['Genome'] = complete_binary['Genome'].str.strip('"')

    new = complete_binary[complete_binary['Genome'].str.match(genome2)]
    distance = new['Distance']
    return distance
