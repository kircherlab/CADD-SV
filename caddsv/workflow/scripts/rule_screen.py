import pandas as pd
import sys
import pysam
from collections import Counter

input_bed = pd.read_table(sys.argv[1], header=None)
input_anno = pysam.TabixFile(sys.argv[2])
contigs = set(input_anno.contigs)

# Define the categories to track
categories_to_track = ["dELS", "pELS", "CA", "CA-CTCF", "CA-H3K4me3", "CA-TF", "TF", "PLS"]

for i, row in input_bed.iterrows():
    chrom, start, end = str(row[0]), row[1], row[2]
    counts = Counter()

    if chrom in contigs:
        for region in input_anno.fetch(chrom, start, end, parser=pysam.asBed()):
            fields = str(region).strip().split("\t")
            if len(fields) >= 6:
                category = fields[5]
                counts[category] += 1

    # Output format: chrom, start, end, each category count, ELS_total, CA_total
    row_counts = [str(counts.get(cat, 0)) for cat in categories_to_track]
    els_total = counts.get("dELS", 0) + counts.get("pELS", 0)
    ca_total = counts.get("CA", 0) + counts.get("CA-CTCF", 0) + counts.get("CA-H3K4me3", 0) + counts.get("CA-TF", 0)
    
    print("\t".join(map(str, [chrom, start, end] + row_counts + [els_total, ca_total])))
