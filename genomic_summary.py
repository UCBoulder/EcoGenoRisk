import openpyxl
import os
import re
import numpy


def genome_extractor(dest):
    EC_num = openpyxl.load_workbook("EC_num")
    sh = EC_num.active
    col = 1
    # name1 = str("row" + i)
    # name2 = str("row" + c)
    dimensions = EC_num.shape
    name1 = numpy.zeros(dimensions)

    for item in os.lisrdir(dest):
        if item.endswith("_matches"):
            gn = open(item, 'r')
            line = gn.readline()
            for row in range(1, sh.max_row + 1):
                if re.search(sh["A" + row].value, line) != "None":
                    print = ("")
                    name1[row, col] = 1
                    break
                numpy.append(name1, numpy.zeros(dimensions), axis=1)
            col += 1
        else:
            break
