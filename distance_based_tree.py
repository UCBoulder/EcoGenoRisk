# This function completes a genome to genome comparitive analysis.
import os
import re
import numpy as np

directory = '/home/anna/PycharmProjects/pythonProject/DIAMOND_matches'
os.chdir(directory)
# opens DIAMOND output, completes difference/similarity
# read_file = open('archaea_big_matrix.csv', 'r')
# matrix_summary = read_file.readlines()
# head = open('archaea_big_matrix.csv', 'r')
# sim_matrix = []
# dif_matrix = []
# genome_sim = []
# genome_diff = ['Genome']
# vertical = [' ']
# horizontal = [' ']
# for genome in head:
#     genome_name = genome.split()[0]
#     if genome_name.startswith("Name_of_Genome"):
#         print(" ")
#     else:
#         horizontal.append(genome_name)
#         vertical = np.vstack([vertical, genome_name])
# for line in matrix_summary:
#     if re.search(line, matrix_summary[0]):
#         print(" ")
#     else:
#         standard = line.split()
#         length = len(standard)
#         genome_vertical = standard[0]
#         for item in matrix_summary:
#             if re.search(item, matrix_summary[0]):
#                 print(" ")
#             else:
#                 comparison = item.split()
#                 sim_sum = 0
#                 diff_sum = 0
#                 counter = 0
#                 similarity = 0
#                 for cell_value in comparison:
#                     if cell_value == '1' or cell_value == '0':
#                         counter += 1
#                         diff_sum = diff_sum + abs(int(cell_value) - int(standard[counter]))
#                 dif_matrix.append(diff_sum)
#         genome_diff = np.vstack([genome_diff, diff_sum])
# print("Let's see the shit show")
# #genome_diff = np.stack([horizontal, genome_diff])
# #genome_diff = vertical.append(genome_diff)
# np.savetxt("Genome_to_genome_comparison", genome_diff, fmt='%s')
horizontal = ['']
comparison = []
vertical=['']
with open('archaea_big_matrix.csv', 'r') as header:
    genome = header.readline()
    for genome in header:
        genome_name = genome.split()[0]
        if genome_name.startswith("Name_of_Genome"):
            print(" ")
        else:
            horizontal.append(genome_name)
            vertiical = np.vstack([vertical, genome_name])

matrix_2 = open('archaea_big_matrix.csv', 'r')
matrix_2_compare = matrix_2.readlines()
matrix_s=open('archaea_big_matrix.csv','r')
matrix_summary = matrix_s.readlines()
valid=0
for line in matrix_summary:
    print("OOWEEEE")    # currently prints 408 of those statements
    split_space = line.split()
    if split_space[0].startswith("Name_of_Genome"):
        print()
        valid=1
    else:
        for line2 in matrix_2_compare:
            compare_line = line2.split()
            diff_sum = 0
            counter = 1
            valid = 0
            for bit in compare_line:
                if re.search("GCF", bit) or re.search(".",bit):
                    valid = 1
                else:
                    diff_sum = diff_sum + abs(int(bit) - int(split_space[counter]))
                    counter += 1
                    print(diff_sum)
            comparison.append(diff_sum)     #after this line runs, the matrix doesn't start a new row and just keeps adding more numbers
    if valid==0:
        horizontal = np.vstack([horizontal, comparison])


np.savetxt("Genome_to_genome_comparison_check", horizontal, fmt='%s')
#if os.path.abspath("Genome_to_genome_comparison_check"):
    #yeet = np.vstack([vertical, horizontal])
    #np.savetxt("Genome_to_genome_comparison_w/_headers", vertical, fmt='%s')
