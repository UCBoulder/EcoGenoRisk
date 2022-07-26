import pandas as pd

pseudo_check = pd.read_csv('/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/Pseudomonas Testing/pseudonomonas_verification.txt', delimiter = '\t')
vertical_pseudo = pseudo_check.T
vertical_pseudo.to_csv('vertical check', sep= '\t')
#
# import subprocess
# print('Sending data to R script')
# # doc_name_1 is sent from function read_in_binary_matrix()
# # doc_name_2 is sent from function tax_clustering()
# # doc_name_3 is sent from function tax_clustering()
# # doc_name_4 is sent from function calculating_distance()
# def to_r_data_processing(doc_name_1, doc_name_2, doc_name_3, doc_name_4, location):
#     #doc_name_1 = 'complete_binary_matrix.txt'
#     #doc_name_2 = 'ec_space_class.txt'
#     #doc_name_3 = '/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/taxonomy.tsv'
#     #doc_name_4 = 'synbio_test_Class_clustered_Distance_Matrix_Combined.txt'
#     #location = '/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/synbio_test/DIAMOND_matches'
#     command = 'Rscript'
#     path2script = '/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/_7_pca_dendro.R'
#     args= [location+'/'+doc_name_1, location+'/'+doc_name_2, doc_name_3, location+'/'+doc_name_4]
#     cmd =[command, path2script] +args
#     #subprocess.check_output(cmd, universal_newlines=True)
#     subprocess.call(['Rscript', '_7_pca_dendro.R', doc_name_1, doc_name_2, doc_name_3, doc_name_4, location])
#     return 'Sent'
#
# to_r_data_processing('complete_binary_matrix.txt','ec_space_class.txt','/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/taxonomy.tsv','synbio_test_Class_clustered_Distance_Matrix_Combined.txt','/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/synbio_test/DIAMOND_matches')
#
# print('Analysis complete.')
