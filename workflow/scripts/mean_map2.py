import pandas as pd
import numpy as np
import sys
import pysam
input_bed = pd.read_table(sys.argv[1], header=None)
input_anno = pysam.TabixFile(sys.argv[2])
contigs = set(input_anno.contigs)
#print(contigs)
for i, row in input_bed.iterrows():
    window = []
    if str(row[0]) in contigs:
        #print("goes on contigs")
        for gc in input_anno.fetch(str(row[0]), row[1], row[2], parser=pysam.asBed()):
            window.append(float(gc[3]))
        if len(window )>0:
            row[3] = np.nanmean(window)
        else:
            row[3] = "."
        print(str(row[0]), row[1], row[2], row[3], sep="\t")
    else:
        print(str(row[0]), row[1], row[2], ".", sep="\t")