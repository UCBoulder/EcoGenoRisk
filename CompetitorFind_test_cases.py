import pandas as pd
import numpy as np
import os as os
import fnmatch
import shutil
import random
import scipy as scipy
import sklearn
from scipy.spatial.distance import pdist, squareform
from sklearn.metrics.pairwise import pairwise_distances


# ======================================================================================================================#

# Creation of similar, different, and random genome pairs for CompetitorFind
# EC Number Binary Matrix is imported to construct the distance matrix which is used to determine similar and different pairs
# Three documents are created that contain the genome pairs
def test_case_genome_isolation():
    # EC BSM is imported
    os.chdir('/home/anna/Desktop/EcoGenoRisk/HazID/CompetitorFind')
    complete_binary = pd.read_csv(
        '/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/combined_binary_summary_matrix_2023_3_21.csv',
        delimiter='\t', header=0, index_col=0)
    # All genome names are extracted
    name_of_all_genomes = complete_binary.index
    # Distance matrix is created from the EC BSM
    distances_parallel = pairwise_distances(X=complete_binary, metric='euclidean', n_jobs=18)
    print('Distance matrix is complete')  #: \, distances_parallel)
    # Similar genomes are found by isolating pairs that have a functional distance greater than zero and less than four
    most_similar_array = distances_parallel[(distances_parallel > 0) & (distances_parallel < 4)]
    # Unique similar genome pairs are created so no pair will be presented twice
    most_sim_uni = np.unique(most_similar_array)
    # Distances are sorted based on their distance, from low to high
    most_sim_uni.sort()
    # Top 10 genomes are selected
    most_similar_array = most_sim_uni[0:10]
    # Different genomes are found by isolating pairs that have a functional distance greater than 18
    most_different_array = distances_parallel[(distances_parallel > 18)]
    # Unique different genome pairs are created so no pair will be presented twice
    most_dif_uni = np.unique(most_different_array)
    # Distances are sorted based on their distance, from low to high
    most_dif_uni.sort()
    # Bottom 10 genomes are selected (largest 10 distances)
    most_different_array = most_dif_uni[-10:]
    # Column names are created for each genome pair dataframe
    genomes_s = ['Similar_Genome_1', 'Similar_Genome_2']
    genomes_d = ['Different_Genome_1', 'Different_Genome_2']
    genome_r = ['Random_Genome_1', 'Random_Genome_2']

    # Random 1000 pairs are chosen
    random_set1 = random.choices(name_of_all_genomes, k=1000)
    random_set2 = random.choices(name_of_all_genomes, k=1000)
    random_stack = np.column_stack((random_set1, random_set2))
    genome_r = np.vstack((genome_r, random_stack))

    # For each distance, corresponding genome names are found for the similar category
    # Once genomes are found, they are horizontally appended to the summary matrix
    for s_value in most_similar_array:
        [i, j] = np.where(distances_parallel == s_value)
        addition = [name_of_all_genomes[i[0]], name_of_all_genomes[j[0]]]
        genomes_s = np.vstack((genomes_s, addition))

    # For each distance, corresponding genome names are found for the similar category
    # Once genomes are found, they are horizontally appended to the summary matrix
    for d_value in most_different_array:
        [k, l] = np.where(distances_parallel == d_value)
        addition = [name_of_all_genomes[k[0]], name_of_all_genomes[l[0]]]
        genomes_d = np.vstack((genomes_d, addition))
    # Distance matrix is converted into a dataframe, column and row names are added
    distances_parallel_df = pd.DataFrame(distances_parallel)
    distances_parallel_df.columns = name_of_all_genomes
    distances_parallel_df.index = name_of_all_genomes

    # Turn on to save the most similar/different genomes list, as well as the overall labelled matrix
    np.savetxt('scoring_tests_similar_genomes.txt', genomes_s, fmt='%s', delimiter='\t')
    np.savetxt('scoring_tests_different_genomes.txt', genomes_d, fmt='%s', delimiter='\t')
    np.savetxt('scoring_tests_random_genomes.txt', genome_r, fmt='%s', delimiter='\t')
    np.savetxt('scoring_tests_distance_matrix_labelled.txt', distances_parallel, fmt='%s', delimiter='\t')
    # Returns the genome list for each category (similar, different, random), complete binary matrix (for identifying
    # each category genomes EC BSM), and the distance matrix (for output file that locates the distance between genomes)
    return genomes_s, genomes_d, genome_r, complete_binary, distances_parallel_df


# ======================================================================================================================#
# Requires the two column dataset that contains the genome pairs, complete binary EC BSM, directory where CompetitorFind
# files are saved, distance matrix
def hijack_competitor_find(genome_set, complete_binary, pathway, type, distances_parallel):
    # Removes the first row for formatting
    genome_set = np.delete(genome_set, 0, axis=0)
    # Iterates through genome list, isolates the rows where genome A and genome B are located
    for genome in genome_set:
        # Genome A and Genome B are isolated. Genome A will be treated as "top match" while Genome B will be used as
        # "synbio". Name designation is arbitrary and results will not change is Genome A is synbio and vice versa
        [genomeA, genomeB] = genome
        # EC BSM is located for genome A
        top_match_bsm = pd.DataFrame(complete_binary.loc[genomeA], columns=[genomeA])
        # Index name is changed in order to merge later downstream
        top_match_bsm.index.name = 'EC-Number'
        # EC BSM is located for genome B
        synbio_bsm = pd.DataFrame(complete_binary.loc[genomeB], columns=[genomeB])
        # Index name is changed in order to merge later downstream
        synbio_bsm.index.name = 'EC-Number'
        # Top match and synbio EC BSMs are merged together so the absence/presence matrix is present for both
        different_ECs = pd.merge(top_match_bsm, synbio_bsm, on='EC-Number', how='outer')
        # Isolates where the EC number absence/presence does not match between the genomes
        different_ECs = different_ECs[different_ECs[genomeA] != different_ECs[genomeB]]
        # Turn on to save the different EC number profiles
        # different_ECs.to_csv(genomeA + '_' + genomeB + '_'+type+'_EC_profile.txt', header=True, index=True, sep='\t')
        # Isolates all EC Numbers that are present for top match for overlap analysis
        top_match_bsm = top_match_bsm[(top_match_bsm.loc[:, genomeA] != 0)]
        # Isolates all EC Numbers that are present for synbio for overlap analysis
        synbio_bsm = synbio_bsm[(synbio_bsm.loc[:, genomeB] != 0)]
        # Sends the genome names, individual EC BSMs, the list of EC numbers that different between genomes,
        # the list of present EC numbers for each organism
        main(genomeA, genomeB, pathway, synbio_bsm, top_match_bsm, different_ECs, type, distances_parallel)
    return "Testing all test cases is complete"


# ======================================================================================================================#

def main(genomeA, genomeB, pathway, synbio_bsm, top_match_bsm, different_ECs, type, distances_parallel):
    # Creates the genome pair specific tag, used for saving documents so no document will be overwritten
    sp_name = genomeA + '_n_' + genomeB + '_' + type
    # Sends the synbio and top match EC BSM to find the total substrate overlap
    individual_genome_rxns = substrate_changes_synbio_v_chassis(sp_name, pathway, synbio_bsm, top_match_bsm)
    # Sends the list of different EC numbers to find the substrate overlap due to pathway modification
    competition = competition_modified_pathway(sp_name, pathway, different_ECs)
    # Finds the duplicitious compounds shared between the total substrate overlap and the modified pathway
    non_repeated_InChIKeys = set(individual_genome_rxns) ^ set(competition)
    non_repeated_InChIKeys = pd.DataFrame(non_repeated_InChIKeys, columns=['InChI-Key'])
    # Calculates the mutualistic compounds where top match reactants and synbio products are considered
    mutualism1 = mutualism1_modified_pathway(sp_name, pathway, synbio_bsm, top_match_bsm)
    # Calculates the mutualistic compounds where top match products and synbio reactants are considered
    mutualism2 = mutualism2_modified_pathway(sp_name, pathway, synbio_bsm, top_match_bsm)
    # Sends the substrate lists to assign a total score
    scoring(individual_genome_rxns, competition, non_repeated_InChIKeys,
            mutualism1, mutualism2, genomeA, genomeB, distances_parallel, sp_name)
    return "Analysis for genome set is complete... "


# ======================================================================================================================#

def scoring(individual, competition, nonrepeated, mutualistic1, mutualistic2, genomeA, genomeB, distances_parallel,
            sp_name):
    # Different weights are applied based on the average number of substrates that are found for each list type
    score = individual.shape[0] / 100 + competition.shape[0] / 10 + nonrepeated.shape[0] - (mutualistic1.shape[0] +
                                                                                            mutualistic2.shape[0]) / 200

    data = np.array([['Individual Reactions', 'Competitive Pathway', 'Unique Compounds', 'Mutualistic1 Pathway',
                      'Mutualistic 2 Pathway', 'Overall Score', 'Distance Score'],
                     [individual.shape[0], competition.shape[0], nonrepeated.shape[0], mutualistic1.shape[0],
                      mutualistic2.shape[0], score, distances_parallel.loc[genomeA, genomeB]]])
    print(data)
    overview = pd.DataFrame(data)
    overview.to_csv(sp_name + '_competitorFind_Analysis.txt', header=True, index=True, sep='\t')


# ======================================================================================================================#
# Moves the files to the category-specific folder. E.I all files that are evaluated for similar genomes are moved to
# one folder
def file_management(pathway, type):
    os.makedirs(type)
    for file in pathway:
        if type in file:
            shutil.move(os.path.join(pathway, file), os.path.join(pathway, type))
    print('Everything has been moved')
    return "File management complete"


# ======================================================================================================================#

print('Test case building is beginning ')
pathway = '/home/anna/Desktop/EcoGenoRisk/HazID/CompetitorFind'

# Initiates test case building
[genome_set1, genome_set2, genome_set_3, complete_binary, distances_parallel] = test_case_genome_isolation()
# Imports all necessary functions
from separatefunctions import substrate_changes_synbio_v_chassis, competition_modified_pathway, \
    mutualism1_modified_pathway, mutualism2_modified_pathway, inchikey_to_conventional_names

# Assigns category types and activates the appropriate functions
type = 'similar'
hijack_competitor_find(genome_set1, complete_binary, pathway, type, distances_parallel)
file_management(pathway, type)
print('Analysis for similar genome is complete, beginning analysis for set of different genomes')
type = 'different'
hijack_competitor_find(genome_set2, complete_binary, pathway, type, distances_parallel)
file_management(pathway, type)
print('Analysis for different genome is complete, beginning analysis for set of random genomes')
type = 'random'
hijack_competitor_find(genome_set_3, complete_binary, pathway, type, distances_parallel)
file_management(pathway, type)
print('Analysis is complete...')
