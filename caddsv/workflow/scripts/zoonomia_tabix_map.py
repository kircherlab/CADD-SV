import sys

import pandas as pd
import pysam


def bedtools_float(value):
    return f"{value:.10g}"


input_bed = pd.read_table(sys.argv[1], header=None)
input_anno = pysam.TabixFile(sys.argv[2])
contigs = set(input_anno.contigs)

for _, row in input_bed.iterrows():
    chrom, start, end = str(row[0]), int(row[1]), int(row[2])
    values = []

    if chrom in contigs and end > start:
        for line in input_anno.fetch(chrom, max(0, start - 1), end):
            fields = line.strip().split("\t")
            anno_start, anno_end = int(fields[1]), int(fields[2])
            if anno_start < end and anno_end > start:
                values.append(float(fields[3]))

    if values:
        print(row[0], row[1], row[2], bedtools_float(max(values)), bedtools_float(min(values)), bedtools_float(sum(values) / len(values)), sep="\t")
    else:
        print(row[0], row[1], row[2], ".", ".", ".", sep="\t")
