# Creates summary matrix of all complete genomes and their EC numbers
import os
import re
import numpy as np


def genome_extractor(diamond, dom):
    os.chdir(diamond)
    ec_open = np.loadtxt("EC_num.csv", dtype='str')  # change to automatic download of EC_num
    print(ec_open)
    big_matrix = ["Name_of_Genome"]
    input_name = input("Save EC binary as (no spaces or capital letters): ")
    file_name= input_name+".csv"
    new_dir = diamond + '/'+file_name
    if os.path.exists(new_dir):
        print("Summary Matrix exists")  # checks if matrix exists
    else:
        for ec_force in ec_open:  # creates a title matrix  (genome name vs. EC number)
            big_matrix.append(ec_force)
        for item in os.listdir(diamond):
            if item.endswith("_matches"):
                genome = [item]  # name of DIAMOND file
                for ec in ec_open:
                    fin = open(item, 'r')  # opens DIAMOND output
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
        np.savetxt(file_name, big_matrix, fmt='%s')
        return new_dir, file_name


def distance_based_tree(ec_space, name_ecfile, dom):
    # This function completes a genome to genome comparative analysis.
    os.chdir(ec_space)
    horizontal1 = []
    vertical = ['']
    with open(name_ecfile, 'r') as header:
       # genome = header.readline()
        for genome in header:
            genome_name = genome.split()[0]
            if genome_name.startswith("Name_of_Genome"):
                print(" ")
            else:
                horizontal1.append(genome_name)  # creates vertical genome labels
                vertical = np.vstack([vertical, genome_name])  # creates horizontal genome labels
    vertical = np.vstack([vertical, ''])
    header.close()
    matrix_2 = open(file_name, 'r')
    matrix_2_compare = matrix_2.readlines()
    matrix_s = open(file_name, 'r')
    matrix_summary = matrix_s.readlines()
    valid = 0
    comparison = []
    for line in matrix_summary:
        print("Processing new genome>>>")
        split_space = line.split()
        if split_space[0].startswith("Name_of_Genome"):
            print()
            valid = 1
        else:
            for line2 in matrix_2_compare:
                compare_line = line2.split()
                diff_sum = 0
                counter = 1
                valid = 0
                for bit in compare_line:

                    if re.search("GCF", bit) or re.search('\.' , bit) or re.search("Name_of_Genome", bit):
                        valid = 1
                    else:
                        diff_sum = diff_sum + abs(int(bit) - int(split_space[counter]))
                        counter += 1
                        valid = 0
                if valid == 0:
                    comparison.append(
                        diff_sum)  # after this line runs, the matrix doesn't start a new row and just keeps adding more numbers
        if valid == 0:
            horizontal = np.vstack([horizontal, comparison])
            comparison = []

    save_mat = np.hstack([vertical, horizontal])
    np.savetxt("Genome_to_genome_difference_comparison", save_mat, fmt='%s')
    header.close()
