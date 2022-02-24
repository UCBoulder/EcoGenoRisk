# User input. Choose synbio file
import subprocess
import os
import shutil
import numpy as np
import re

print("Welcome to EcoGeno! Enter your pathway for the synbio FASTA file: ")
sb_path = input()
print("Enter where you want results saved: ")
desired_location = input()
finding_name = sb_path.split("/")
group_num = len(finding_name)
sb_name = finding_name[group_num - 1]
matches = sb_name + "_matches"
new = desired_location + '/' + matches

# DIAMOND run --complete
if os.path.abspath(new):
    print(" ")
else:
    blastp = ['diamond', 'blastp', '-d', '/home/anna/PycharmProjects/pythonProject/DIAMOND_matches/reference.dmnd',
              '-q',
              sb_path,
              '-o', matches, '--max-target-seqs', '1', '--outfmt', '6']  # completes DIAMOND search by using full path
    subprocess.run(blastp)
    shutil.move(os.path.abspath(matches), desired_location)

# Summary Matrix Rendering--complete
directory = '/home/anna/PycharmProjects/pythonProject/DIAMOND_matches'
os.chdir(directory)
ec_open = np.loadtxt("EC_num.csv", dtype='str')
big_matrix = ['Name_of_genome']
name = sb_name + "_big_matrix.csv"
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
np.savetxt(name, big_matrix, fmt='%s')

# Matrix comparison--in progress
# Method 1: Difference Based Sum --complete

print("Domain of synbio? ")  # finding domain and corresponding matrix
domain = input()
lc_domain = domain.lower()
if lc_domain == 'archaea':
    archaea = open('archaea_big_matrix.csv', 'r')
    summary_mat = archaea.readlines()
if lc_domain == 'bacteria':
    bacteria = open('bacteria_big_matrix', 'r')
    summary_mat = bacteria.readlines()
if lc_domain == 'fungi':
    fungi = open('fungi_big_matrix', 'r')
    summary_mat = fungi.readlines()
else:
    other = open(name, 'r')
sb_matrix = open(name, 'r')
content = sb_matrix.readlines()
sm = content[1].split()  # splits second line by tab
m1_comparison = ["Genome", "Difference_Based_Sum"]
m1_overall = ['-', '-']

for line in summary_mat:
    counter = 0
    diff_sum = 0
    if re.search(line, summary_mat[0]):
        print("NA")
    else:
        bm_split = line.split()
        for z in bm_split:
            if z == '0' or z == '1':
                counter += 1
                diff_sum = diff_sum + abs(int(sm[counter]) - int(z))
            m1_overall = [bm_split[0], diff_sum]
        m1_comparison = np.vstack((m1_comparison, m1_overall))
save_name = sb_name + "_diffsum_comparison.csv"
np.savetxt(save_name, m1_comparison, fmt='%s')

# Method 2: Similarity Based --in progress
m2_comparison = ["Genome", "Difference_Based_Sum"]
m2_overall = ['-', '-']
for line in summary_mat:
    count = 0
    sim_sum = 0
    similarity = 0
    if re.search(line, summary_mat[0]):
        print("NA")
    else:
        bm_split = line.split()
        for z in bm_split:
            if z == '0' or z == '1':  # if big matrix value is a digit
                count += 1
                if z == '1' and sm[count] == '1':
                    similarity = 2
                    sim_sum = sim_sum + similarity
                    print(sim_sum)
                if z == '0' and sm[count] == '0':
                    similarity = 0
                    sim_sum = sim_sum + similarity
                else:
                    similarity = -1
                    sim_sum = sim_sum + similarity
            overall = [bm_split[0], sim_sum]
        m2_comparison = np.vstack((m2_comparison, overall))
m2_save_name = sb_name + "_simsum_comparison.csv"
np.savetxt(m2_save_name, m2_comparison, fmt='%s')

# Moving all files
shutil.move(os.path.abspath(name), desired_location)
shutil.move(os.path.abspath(save_name), desired_location)
shutil.move(os.path.abspath(m2_save_name), desired_location)
print("That's all folks!")
