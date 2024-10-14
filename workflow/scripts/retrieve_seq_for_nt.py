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
        upstr = hg38fasta.fetch(row[0], row[1]-2400, row[1])
        downstr = hg38fasta.fetch(row[0], row[2], row[2]+2400)
        # Check the length of row[4]
        if len(row[4]) > 25200:
            # If the length of row[4] is greater than 25200 characters
            first_part = row[4][:12450]  # First 12450 characters
            last_part = row[4][-12450:]  # Last 12450 characters
            row4_modified = first_part + 'N' * 50 + last_part  # Insert 50 'N's in between
        else:
            # If the length of row[4] is less than or equal to 25200 characters
            row4_modified = row[4]
        seq = upstr + row4_modified + downstr
        print(row[0], row[1], row[2], row[3], seq.upper(), sep="\t")
    elif (row[3] == "DEL"):
        upstr = hg38fasta.fetch(row[0], row[1]-2400, row[1])
        downstr = hg38fasta.fetch(row[0], row[2], row[2]+2400)
        seq = upstr + downstr
        print(row[0], row[1], row[2], row[3], seq.upper(), sep="\t")
    elif (row[3] == "DUP"):
        row[4] = hg38fasta.fetch(row[0], row[1], row[2])
        upstr = hg38fasta.fetch(row[0], row[2]-2400, row[2])
        downstr = hg38fasta.fetch(row[0], row[2], row[2]+2400)
        if len(row[4]) > 25200:
            # If the length of row[4] is greater than 25200 characters
            first_part = row[4][:12450]  # First 12450 characters
            last_part = row[4][-12450:]  # Last 12450 characters
            row4_modified = first_part + 'N' * 50 + last_part  # Insert 50 'N's in between
        else:
            # If the length of row[4] is less than or equal to 25200 characters
            row4_modified = row[4]
        seq = upstr + row4_modified + downstr
        print(row[0], row[1], row[2], row[3], seq.upper(), sep="\t")
'''        
for i, row in input.iterrows():
    upstr = hg38fasta.fetch(row[0], row[1]-2400, row[1])
    downstr = hg38fasta.fetch(row[0], row[2], row[2]+2400)
    seq = upstr + downstr
    print(row[0], row[1], row[2], seq.upper(), sep="\t")
'''