import openpyxl
import os
import re
import numpy
from numpy import savetxt


def genome_extractor(dest):
    os.chdir("/home/anna/PycharmProjects/pythonProject/DIAMOND_matches")
    ec_num = openpyxl.load_workbook("EC_num.xlsx")
    sh = ec_num.active
    col = 1
    row_num = sh.max_row
    # dimensions = ec_num.shape
    # (row_num, colm) = dimensions
    name1 = numpy.zeros((row_num, 1))
    # f = open('genome_output' , 'w')

    for item in os.listdir(dest):
        if item.endswith("_matches"):
            print(item)
            gn = open(item, 'r')
            line = gn.readline()
            for row in range(1, row_num):
                if re.search(sh["A" + str(row)].value, line) != "None":
                    name1[row, col] = 1
                    # numpy.savetext()
            col += 1
        name1 = numpy.append(name1, numpy.zeros((row_num, 1)), axis=1)
    print(line)
    savetxt('data.cvs',name1, delimiter=',')
