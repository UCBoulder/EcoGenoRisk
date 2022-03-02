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
horizontal = []
comparison = []
with open('archaea_big_matrix.csv', 'r') as header:
    genome = header.readline()
    for genome in header:
        genome_name = genome.split()[0]
        if genome_name.startswith("Name_of_Genome"):
            print(" ")
        else:
            horizontal.append(genome_name)
            if genome_name=="GCF_000968395.2_ASM96839v2_protein_matches":
                vertical=["GCF_000968395.2_ASM96839v2_protein_matches"]
            else:
                vertical = np.vstack([vertical, genome_name])
print(horizontal)
print(vertical)
matrix_2 = open('archaea_big_matrix.csv', 'r')
matrix_2_compare = matrix_2.readlines()
with open('archaea_big_matrix.csv', 'r') as matrix_summary:
    line = matrix_summary.readline()
    print("OOWEEEE")
    counter = 1
    for line in matrix_summary:
        split_space = line.split()
        for line2 in matrix_2_compare:
            compare_line = line2.split()
            diff_sum = 0
            if compare_line == '1' or compare_line == '0':
                for i in compare_line:
                    if i == '1' or i == '0':
                        diff_sum = diff_sum + abs(int(i) - int(split_space[counter]))
                comparison.append(diff_sum)
            print(diff_sum)
        counter += 1
    horizontal = np.vstack([horizontal, comparison])
np.savetxt("Genome_to_genome_comparison_check", horizontal, fmt='%s')
if os.path.abspath("Genome_to_genome_comparison_check"):
    yeet = np.vstack([vertical, horizontal])
    np.savetxt("Genome_to_genome_comparison_w/_headers", vertical, fmt='%s')
