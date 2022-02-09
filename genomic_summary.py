import numpy as np
import openpyxl
import os
import re
import numpy as np
from numpy import savetxt


def genome_extractor(dest):
    os.chdir("/home/anna/PycharmProjects/pythonProject/DIAMOND_matches")
    ec_open = np.loadtxt("EC_num.csv", dtype='str')
    print(ec_open)
    name1 = np.zeros(len(ec_open))
    big_matrix = ["Name of Genome"]
    for ec_force in ec_open:    #creates a title matrix  (genome name vs. EC number)
        big_matrix.append(ec_force)
    print(big_matrix)
    print(name1)
    counter = 0
    for item in os.listdir(dest):
        #if counter < 2: #allows only two genomes to be processed
        #   counter += 1
        if item.endswith("_matches"):
            genome = [item]
            for ec in ec_open:
                fin = open(item, 'r')
                ec_stand_in = np.zeros(len(ec_open))
                ec_now = 0
                for line in fin:
                    x = line.split('\t')
                    #print(x[1])
                    y = x[1].split("?")
                    #print(y[1])
                    z = y[1].split(";_")
                    additional = len(z)
                    #print(additional)
                    #print(z)

                    if re.fullmatch(ec, y[1]) is not None:  # currently looks if the numbers are present and not for a direct match
                        ec_now = 1
                    # Extend list by that EC
                    for i in z:
                        if re.fullmatch(ec, i) is not None:
                            ec_now = 1

                genome.append(ec_now)
            big_matrix = np.vstack([big_matrix, genome])
    print(big_matrix)
    np.savetxt('big_matrix.cvs', big_matrix, fmt='%s')
