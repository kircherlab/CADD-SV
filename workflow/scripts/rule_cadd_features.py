import pandas as pd
import numpy as np
import sys
import pysam

input_bed = pd.read_table(sys.argv[1], header=None)
input_anno = pysam.TabixFile(sys.argv[2])
contigs = set(input_anno.contigs)

positive_only_indices = [3, 4, 13, 14]
may_have_negatives_indices = [i for i in range(3, 15) if i not in positive_only_indices]

for _, row in input_bed.iterrows():
    chrom, start, end = row[0], row[1], row[2]
    values = {i: [] for i in range(3, 15)}

    if chrom in contigs:
        for line in input_anno.fetch(str(chrom), start, end):
            fields = line.strip().split('\t')
            for i in range(3, 15):
                if i < len(fields):
                    val = 0.0 if fields[i] == "." else float(fields[i])
                    values[i].append(val)

        if not any(values[i] for i in values):
            print(chrom, start, end, *["."] * 32, sep="\t")
        else:
            output = [chrom, start, end]
            for i in range(3, 15):
                col = values[i]
                if not col:
                    stats = ["."] * (2 if i in positive_only_indices else 3)
                elif i in positive_only_indices:
                    stats = [max(col), sum(col)]
                else:
                    stats = [min(col), max(col), np.sum(np.abs(col))]
                output.extend(stats)
            print(*output, sep="\t")
    else:
        print(chrom, start, end, *["."] * 32, sep="\t")
