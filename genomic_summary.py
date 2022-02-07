import numpy as np
import openpyxl
import os
import re
import numpy as np
from numpy import savetxt


def genome_extractor(dest):
    os.chdir("/home/anna/PycharmProjects/pythonProject/DIAMOND_matches")
    ec_open = np.loadtxt("EC_num.csv",dtype='str')
    #ec_table = ec_open.readlines()
    #print(len(ec_table))
    col = 1
    #row_num = len(ec_table)
    print(ec_open)
    # dimensions = ec_num.shape
    # (row_num, colm) = dimensions
    name1 = np.zeros(len(ec_open))

    big_matrix = ["Name of Genome"]
    for ec_force in ec_open:
        big_matrix.append(ec_force)
    print(big_matrix)

    print(name1)

    counter = 0
    # f = open('genome_output' , 'w')
    for item in os.listdir(dest):
        if counter < 2:
            counter += 1
            if item.endswith("_matches"):
        #    print(item)
        #    gn = open(item, 'r')
        #    line = gn.readline()
        #    print(line)
        #    temp_ec_match = numpy.zeros(row_num)

        #    for row in range(1, row_num):
        #        ec_entry = ec_table[row]
        #        #print(ec_entry.strip())
        #        if re.search(ec_entry.strip(), line):
        #            temp_ec_match[row] = 1
        #            print("match")
        #        else:
        #            temp_ec_match[row] = 0
        #    col += 1
        #    print(temp_ec_match)
        #name1 = numpy.append(name1, temp_ec_match, axis=0)
        #print(line)
                genome = [item]
                for ec in ec_open:
                    fin = open(item,'r')
                    ec_stand_in = np.zeros(len(ec_open))
                    ec_now = 0
                    for line in fin:
                        if re.search(ec,line) is not None:
                            ec_now = 1
                        #Extend list by that EC
                    genome.append(ec_now)
                big_matrix = np.vstack([big_matrix, genome])
        print(big_matrix)
    np.savetxt('big_matrix.cvs',big_matrix, fmt='%s')
