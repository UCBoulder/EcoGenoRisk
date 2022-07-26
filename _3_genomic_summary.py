# Creates summary matrix of all complete genomes and their EC numbers
import os;
import re;
import numpy as np

# This function creates the EC binary summary matrix, by going through all of the DIAMOND search outputs and denoting which
# EC numbers are present in each genome, by using a 1 or 0 for EC number status
# This function requires the location of the DIAMOND (Buchfink et al., 2021) search results found in diamond_impl, name of the domain
# Input requires the list of EC numbers, DIAMOND search results, and domain name
# Output is a space separated textfile that should be num genomes x num of EC num as well as the location of the file
def genome_extractor(diamond, name):
    # Changes the directory to the location of the DIAMOND search outputs
    os.chdir(diamond)
    print(os.getcwd())
    # Opens the list of of EC numbers
    ec_open = np.loadtxt('/home/anna/Desktop/EcoGenoRisk/HazID/NicheOverlap/EC_num.csv',
                         dtype='str')  # change to automatic download of EC_num
    big_matrix = ["Name_of_Genome"]
    # Asks user to input name for EC matrix
    # input_name = input("Save EC binary summary matrix as (no spaces or capital letters): ")
    # Specifies document to be a csv type
    # file_name= input_name+".csv"
    # file_name = "synbio1_big_matrix.csv"
    if name == "":
        file_name = os.path.abspath(diamond).rsplit('/', 9)[7] + '_binary_matrix.txt'
        print(file_name)
    else:
        file_name = name
    new_dir = diamond + '/' + file_name
    # Checks to see if the document already exists using full pathway name
    if os.path.exists(new_dir):
        print("Summary Matrix exists")
        return [new_dir, file_name]
    else:
        for ec_force in ec_open:
            # Creates a horizontal header of all of the EC names
            big_matrix.append(ec_force)
        # Goes through all of the DIAMOND outputs in the folder
        # Goes through all of the output files, each one is opened and read one line at a time. The lines are split to
        # extract the EC numbers found in each line. If the EC number found in the DIAMOND output matches an
        # EC entry in the list, the status is changed from a zero to a one. The binary status is catalogued horizontally
        # for each genome, and following genomes are vertically stacked
        for item in os.listdir(diamond):
            if item.endswith("_matches.tsv"):
                print(item)
                # Finds the name of the DIAMOND output file
                genome = [item]
                print(genome)
                # Iterates through all of the EC numbers (1:8197)
                for ec in ec_open:
                    # Opens individual DIAMOND output files
                    fin = open(item, 'r')
                    # Sets default for EC status is zero, meaning absent
                    ec_now = 0
                    # Takes the first line in the DIAMOND output file and splits it based on tab separation
                    # Takes the second column of the split line, which has EC numbers separated by a ?, :_
                    # Strings splits have a new name assigned to them
                    for line in fin:
                        no_tab = line.split('\t')
                        first_ec = no_tab[1].split("?")
                        separate_ec = first_ec[1].split(";_")
                        # Checks for a full match between the EC number listed in the DIAMOND output and the EC number
                        # found in the separate document
                        if re.fullmatch(ec, first_ec[1]) is not None:  # looks for full match of first EC number
                            ec_now = 1
                        # In the case that there are more than one EC separated by ;, the function iterates through the list
                        # and sees if there is a full match between the listed EC and the list
                        for i in separate_ec:
                            if re.fullmatch(ec, i) is not None:  # looks for full match of any other ECs listed
                                ec_now = 1
                    # 1 or 0 will be appended to the summary matrix for each EC value in the list
                    genome.append(ec_now)
                # Vertical stacking occurs for each genome in the DIAMOND output folder
                big_matrix = np.vstack([big_matrix, genome])
        print(big_matrix)
        # Saves matrix as a text file for further analysis
        np.savetxt(file_name, big_matrix, fmt='%s')
        # Returns the location of the summary matrix and the name of the file
        print(new_dir)
        return [new_dir, file_name]
##=============================================Citations==============================================================##
# Buchfink B, Reuter K, Drost HG, "Sensitive protein alignments at tree-of-life scale using DIAMOND", Nature Methods 18, 366–368 (2021). doi:10.1038/s41592-021-01101-x
