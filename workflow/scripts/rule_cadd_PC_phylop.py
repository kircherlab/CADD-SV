import pandas as pd
import numpy as np
import sys
import pysam
input_bed = pd.read_table(sys.argv[1], header=None)
input_anno = pysam.TabixFile(sys.argv[2])
contigs = set(input_anno.contigs)

for i, row in input_bed.iterrows():
    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []
    if row[0] in contigs:
        for gc in input_anno.fetch(str(row[0]), row[1], row[2], parser=pysam.asBed()):
            col1.append(float(gc[3]))
            col2.append(float(gc[4]))
            col3.append(float(gc[5]))
            col4.append(float(gc[6]))
            col5.append(float(gc[7]))
        if col1 == []:
            print(str(row[0]), row[1], row[2], ".", ".", ".", ".", ".", ".",".", ".", ".", ".", sep="\t")
        else:
            print(str(row[0]), row[1], row[2], max(col1), sum(col1), max(col2), sum(col2), max(col3), sum(col3), max(col4), sum(col4), max(col5), sum(col5), sep="\t")
    else:
        print(str(row[0]), row[1], row[2], ".", ".", ".", ".", ".", ".",".", ".", ".", ".", sep="\t")