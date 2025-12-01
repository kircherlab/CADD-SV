import pandas as pd
import numpy as np
import sys
import pysam

input_bed = pd.read_table(sys.argv[1], header=None)
input_anno = pysam.TabixFile(sys.argv[2])
contigs = set(input_anno.contigs)

for i, row in input_bed.iterrows():
    count = 0
    value_sum = 0
    if str(row[0]) in contigs:
        for region in input_anno.fetch(str(row[0]), row[1], row[2], parser=pysam.asBed()):
            count += 1
            try:
                value_sum += float(region.name)  # region.name refers to the 4th column
            except (ValueError, TypeError):
                continue  # Skip if it's not a number
        print(str(row[0]), row[1], row[2], count, value_sum, sep="\t")
    else:
        print(str(row[0]), row[1], row[2], ".", ".", sep="\t")
