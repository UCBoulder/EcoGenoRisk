# User input. Choose synbio file, for single synbio file
import subprocess
import os
import shutil
import numpy as np
import re
import pandas as pd
from distance_matrix import pass_to_distance

print("Welcome to EcoGeno! Enter the name of your FASTA file: ")
sb_filename = input()
print("Enter where you want results saved: ")
desired_location = input()
print("Enter organism ID, no punctuation: ")
sb_name = input()
matches = sb_filename+"_matches"

# DIAMOND run --complete
if os.path.exists(matches):
    print(" ")
else:
    blastp = ['diamond', 'blastp', '-d', '/home/anna/PycharmProjects/pythonProject/Bacteria/DIAMOND_matches/reference.dmnd',
              '-q',
              desired_location+"/"+sb_filename,
              '-o', matches, '--max-target-seqs', '1', '--outfmt', '6']  # completes DIAMOND search by using full path
    subprocess.run(blastp)
print("DIAMOND Results Completed")
## Summary Matrix Rendering-- Completes the binary matrix using the diamond hits outout (matches)
#

os.chdir(desired_location)
name = sb_name + "_big_matrix.txt"
if os.path.exists(name):
    print(" ")
else:
    ec_open = np.loadtxt("EC_num.csv", dtype='str')  # change to automatic download of EC_num
    big_matrix = ["Name_of_Genome"]
    genome = [sb_name]
    for ec_force in ec_open:  # creates a title matrix  (genome name vs. EC number)
        big_matrix.append(ec_force)
    for ec in ec_open:
        fin = open(matches, 'r')  # opens DIAMOND output
        ec_now = 0
        for line in fin:
            no_tab = line.split('\t')
            first_ec = no_tab[1].split("?")
            separate_ec = first_ec[1].split(";_")
            if re.fullmatch(ec, first_ec[1]) is not None:  # looks for full match of first EC number
                ec_now = 1
            # Extend list by that EC
            for i in separate_ec:
                if re.fullmatch(ec, i) is not None:  # looks for full match of any other ECs listed
                    ec_now = 1
        genome.append(ec_now)
    big_matrix = np.vstack([big_matrix, genome])
    print(big_matrix)
    np.savetxt(name, big_matrix, fmt='%s')

print("EC Binary Scoring is Complete")
distance_list_for_synbio = pass_to_distance(name, sb_name, desired_location)          #Passes the filename of the binary matrix and the ID of organism
print(distance_list_for_synbio)             #Prints distance matrix for the synbio genome
distance_list_for_synbio.to_csv("Distance_List_for_UploadedSynbio_"+sb_name+".txt", sep="\t")
#np.savetxt("DistanceListforSynbio.txt", distance_list_for_synbio, fmt='%s')   #saves file of the synbio distance matrix






# os.chdir(desired_location)
# name = sb_name + "_big_matrix.txt"
# if os.path.exists(sb_name+"_big_matrix.txt"):      #checks if the binary matrix already exists, if so, prints space if not completes the analysis
#     print(" ")
# else:
#     directory = '/home/anna/PycharmProjects/pythonProject/Bacteria/DIAMOND_matches'
#     os.chdir(directory)
#     ec_open = np.loadtxt("EC_num.csv", dtype='str')
#     big_matrix = ['Name_of_genome']
#
#     genome = [sb_name]
#
#     for ec_force in ec_open:  # creates a title matrix  (genome name vs. EC number)
#         big_matrix.append(ec_force)
#     for ec in ec_open:
#         fin = open(matches, 'r')  # opens DIAMOND output
#         ec_now = 0
#         for line in fin:
#             no_tab = line.split('\t')
#             first_ec = no_tab[1].split("?")
#             separate_ec = first_ec[1].split(";_")
#             if re.fullmatch(ec, first_ec[1]) is not None:  # looks for full match of first EC number
#                 ec_now = 1
#             # Extend list by that EC
#             for i in separate_ec:
#                 if re.fullmatch(ec, i) is not None:  # looks for full match of any other ECs listed
#                     ec_now = 1
#         genome.append(ec_now)
#     big_matrix = np.vstack([big_matrix, genome])
#     np.savetxt(name, big_matrix, fmt='%s')          #Note, I converted this output file to a tab separated one
#     #shutil.move(name, desired_location)