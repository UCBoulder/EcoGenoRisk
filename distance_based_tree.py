# # This function completes a genome to genome comparative analysis.
# import os
# import re
# import numpy as np
#
# directory = '/home/anna/PycharmProjects/pythonProject/DIAMOND_matches'
# os.chdir(directory)
# horizontal = []
# vertical=['']
# with open('archaea_big_matrix.csv', 'r') as header:
#     genome = header.readline()
#     for genome in header:
#         genome_name = genome.split()[0]
#         if genome_name.startswith("Name_of_Genome"):
#             print(" ")
#         else:
#             horizontal.append(genome_name)
#             vertical = np.vstack([vertical, genome_name])
# vertical=np.vstack([vertical,''])
# matrix_2 = open('archaea_big_matrix.csv', 'r')
# matrix_2_compare = matrix_2.readlines()
# matrix_s=open('archaea_big_matrix.csv','r')
# matrix_summary = matrix_s.readlines()
# valid=0
# comparison = []
#
# for line in matrix_summary:
#     print("Processing new genome>>>")
#     split_space = line.split()
#     if split_space[0].startswith("Name_of_Genome"):
#         print()
#         valid = 1
#     else:
#         for line2 in matrix_2_compare:
#             compare_line = line2.split()
#             diff_sum = 0
#             counter = 1
#             valid = 0
#             for bit in compare_line:
#
#                 if re.search("GCF", bit) or re.search("\.",bit) or re.search("Name_of_Genome",bit):
#                     valid = 1
#                 else:
#                     diff_sum = diff_sum + abs(int(bit) - int(split_space[counter]))
#                     counter += 1
#                     valid = 0
#             if valid == 0:
#                 comparison.append(diff_sum)     #after this line runs, the matrix doesn't start a new row and just keeps adding more numbers
#     if valid == 0:
#         horizontal = np.vstack([horizontal, comparison])
#         comparison = []
#
# save_mat=np.hstack([vertical,horizontal])
# np.savetxt("Genome_to_genome_difference_comparison", save_mat, fmt='%s')
# header.close()

# This function completes a genome to genome comparative analysis.
import os
import re
import numpy as np

directory = '/home/anna/PycharmProjects/pythonProject/DIAMOND_matches'
os.chdir(directory)
horizontal1 = []
vertical1 = ['']

with open('archaea_big_matrix.csv', 'r') as header:
    genome = header.readline()
    for genome in header:
        genome_name = genome.split()[0]
        if genome_name.startswith("Name_of_Genome"):
            print(" ")
        else:
            horizontal1.append(genome_name)
            vertical1 = np.vstack([vertical1, genome_name])
vertical1 = np.vstack([vertical1, ''])
header.close()
matrix_2 = open('archaea_big_matrix.csv', 'r')
matrix_2_compare = matrix_2.readlines()
matrix_s = open('archaea_big_matrix.csv','r')
matrix_summary = matrix_s.readlines()
valid = 0
comparison = []
matrix=['']

for line in matrix_summary:     # iterates through the lines of main matrix
    print("Processing new genome>>>")
    split_space = line.split()
    if split_space[0].startswith("Name_of_Genome"):
        print('')
        valid = 1
    else:
        for line2 in matrix_2_compare:      # iterates through the lines of matrix to compare
            compare_line = line2.split()
            similarity = 0
            counter = 1
            sim_sum = 0
            for bit in compare_line:
                if re.search("GCF", bit) or re.search("\.", bit) or re.search("Name_of_Genome", bit):
                    valid = 1
                else:
                    valid = 0
                    if bit == '1' and split_space[counter] == '1':
                        similarity = 2
                    elif bit == '0' and split_space[counter] == '0':
                        similarity = 0
                    else:
                        similarity = -1
                    # add = int(bit) + int(split_space[counter])
                    # if add == 1:
                    #     add = add * -1
                    counter += 1
                    sim_sum = similarity + sim_sum
            if valid == 0:
                comparison.append(sim_sum)  # after this line runs, the matrix doesn't start a new row and just keeps adding more numbers
    if valid == 0:
        horizontal1 = np.vstack([horizontal1, comparison])
        comparison = []

np.savetxt("Yeet",horizontal1, fmt='%s')
if len(matrix) != 1:
    matrix=np.vstack([horizontal1, matrix])
    save_mat = np.hstack([vertical1, matrix])
    np.savetxt("Genome_to_genome_similarity_comparison", save_mat, fmt='%s')

