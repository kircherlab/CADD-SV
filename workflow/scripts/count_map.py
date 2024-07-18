import pandas as pd
import numpy as np
import sys
import pysam
input_bed = pd.read_table(sys.argv[1], header=None)
input_anno = pysam.TabixFile(sys.argv[2])
contigs = set(input_anno.contigs)

for i, row in input_bed.iterrows():
    c = 0
    if str(row[0]) in contigs:
        for regions in input_anno.fetch(str(row[0]), row[1], row[2], parser=pysam.asBed()):
            c += 1
        print(str(row[0]), row[1], row[2], c, sep="\t")
    else:
        print(str(row[0]), row[1], row[2], ".", sep="\t")