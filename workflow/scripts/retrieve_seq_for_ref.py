import pandas as pd
import pysam
import sys

def transform_genomic_sequence(sequence):
    transformed_sequence = []
    i = 0
    
    while i < len(sequence):
        if sequence[i] == 'N':
            start = i
            while i < len(sequence) and sequence[i] == 'N':
                i += 1
            n_count = i - start
            transformed_sequence.append('N' * (n_count // 6))
            transformed_sequence.append('A' * (n_count % 6))
        else:
            transformed_sequence.append(sequence[i])
            i += 1
    return ''.join(transformed_sequence)

def trim_span(span_seq):
    if len(span_seq) > 1008:
        return span_seq[:498] + 'N' * 12 + span_seq[-498:]
    return span_seq

# Input
input_file = sys.argv[1]
reference_genome = sys.argv[2]

hg38fasta = pysam.FastaFile(reference_genome)
inputfile = pd.read_table(input_file, header=None)

# Ensure chromosomes have "chr" prefix
bed_inputfile = []
for i, row in inputfile.iterrows():
    if ("chr" not in str(row[0])):
        row[0] = "chr" + str(row[0])
    bed_inputfile.append(row)
bed_inputfile = pd.DataFrame(bed_inputfile)

chrom_list = ["chr" + str(i) for i in list(range(1, 23))] + ["chrX", "chrY"]
input = bed_inputfile[bed_inputfile[0].isin(chrom_list)]

for _, row in input.iterrows():
    chrom, start, end, svtype = row[0], int(row[1]), int(row[2]), row[3]
    upstream_length = min(start, 96)

    if svtype == "INS":
        upstr = hg38fasta.fetch(chrom, start - upstream_length, start)
        downstr = hg38fasta.fetch(chrom, end, end + 96)
        seq = upstr + downstr
    elif svtype in ["DEL", "INV", "DUP"]:
        upstr = hg38fasta.fetch(chrom, start - upstream_length, start)
        span = hg38fasta.fetch(chrom, start, end)
        span = trim_span(span[::-1] if svtype == "INV" else span)
        downstr = hg38fasta.fetch(chrom, end, end + 96)
        seq = upstr + span + downstr
    else:
        continue  # Skip unsupported SV types

    if "N" in seq:
        seq = transform_genomic_sequence(seq)
    print(chrom, start, end, svtype, seq.upper(), sep="\t")
