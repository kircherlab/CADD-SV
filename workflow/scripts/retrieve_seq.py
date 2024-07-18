import pandas as pd
import pysam
import sys

input_file = sys.argv[1]
reference_genome = sys.argv[2]

hg38fasta = pysam.FastaFile(reference_genome)

inputfile = pd.read_table(input_file, header=None)

bed_inputfile = []
for i, row in inputfile.iterrows():
    if ("chr" not in str(row[0])):
        row[0] = "chr" + str(row[0])
    bed_inputfile.append(row)
bed_inputfile = pd.DataFrame(bed_inputfile)

chrom_list = ["chr1", "chr2", "chr3",
              "chr4", "chr5", "chr6",
              "chr7", "chr8", "chr9",
              "chr10", "chr11", "chr12",
              "chr13", "chr14", "chr15",
              "chr16", "chr17", "chr18",
              "chr19", "chr20", "chr21",
              "chr22", "chrX", "chrY"]

input = bed_inputfile[bed_inputfile[0].isin(chrom_list)]


# output = []
for i, row in input.iterrows():
    if (row[3] == "INS"):
        # output.append(row)
        print(row[0], row[1], row[2], row[3], row[4], sep="\t")
    else:
        row[4] = hg38fasta.fetch(row[0], row[1], row[2])
        # output.append(row)
        print(row[0], row[1], row[2], row[3], row[4], sep="\t")
