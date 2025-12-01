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
            
            # Append multiples of 6 as 'N'
            transformed_sequence.append('N' * (n_count // 6))
            # Append remaining as 'A'
            transformed_sequence.append('A' * (n_count % 6))
        else:
            transformed_sequence.append(sequence[i])
            i += 1
    
    return ''.join(transformed_sequence)

def truncate_to_200_tokens(sequence):
    i = 0
    token_count = 0
    output = []

    while i < len(sequence) and token_count < 200:
        if sequence[i] == 'N':
            if token_count + 1 > 200:
                break
            output.append('N')
            token_count += 1
            i += 1
        else:
            # Check if a 6-mer can be formed before the next 'N'
            next_N = sequence.find('N', i)
            end = next_N if next_N != -1 else len(sequence)
            remaining = end - i

            if remaining >= 6:
                if token_count + 1 > 200:
                    break
                output.append(sequence[i:i+6])
                token_count += 1
                i += 6
            else:
                # Add each remaining base as 1 token
                for j in range(remaining):
                    if token_count + 1 > 200:
                        break
                    output.append(sequence[i])
                    token_count += 1
                    i += 1
    return ''.join(output)
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
    if row[1] > 96:
        upstream_length = 96
    else:
        upstream_length = row[1]
    if (row[3] == "INS"):
        upstr = hg38fasta.fetch(row[0], row[1]-upstream_length, row[1])
        downstr = hg38fasta.fetch(row[0], row[2], row[2]+96)
        if pd.isna(row[4]) or not row[4]:
            row[4] = 'NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN'
        # Check the length of row[4]
        if len(row[4]) > 1008:
            # If the length of row[4] is greater than 25200 characters
            first_part = row[4][:498]  # First 12450 characters
            last_part = row[4][-498:]  # Last 12450 characters
            row4_modified = first_part + 'N' * 12 + last_part  # Insert 50 'N's in between
        else:
            # If the length of row[4] is less than or equal to 25200 characters
            row4_modified = row[4]
        seq = upstr + row4_modified + downstr
        if "N" in seq:
            seq = transform_genomic_sequence(seq)
        print(row[0], row[1], row[2], row[3], seq.upper(), sep="\t")
    elif (row[3] == "DEL"):
        upstr = hg38fasta.fetch(row[0], row[1]-upstream_length, row[1])
        downstr = hg38fasta.fetch(row[0], row[2], row[2]+96)
        seq = upstr + downstr
        if "N" in seq:
            seq = transform_genomic_sequence(seq)
        print(row[0], row[1], row[2], row[3], seq.upper(), sep="\t")
    elif (row[3] == "INV"):
        row[4] = hg38fasta.fetch(row[0], row[1], row[2])
        upstr = hg38fasta.fetch(row[0], row[1]-upstream_length, row[1])
        downstr = hg38fasta.fetch(row[0], row[2], row[2]+96)
        if len(row[4]) > 1008:
            # If the length of row[4] is greater than 25200 characters
            first_part = row[4][:498]  # First 12450 characters
            last_part = row[4][-498:]  # Last 12450 characters
            row4_modified = first_part + 'N' * 12 + last_part  # Insert 50 'N's in between
        else:
            # If the length of row[4] is less than or equal to 25200 characters
            row4_modified = row[4]
        seq = upstr + row4_modified[::-1] + downstr
        if "N" in seq:
            seq = transform_genomic_sequence(seq)
        print(row[0], row[1], row[2], row[3], seq.upper(), sep="\t")
    elif (row[3] == "DUP"):
        row[4] = hg38fasta.fetch(row[0], row[1], row[2])
        upstr = hg38fasta.fetch(row[0], row[2]-96, row[2])
        downstr = hg38fasta.fetch(row[0], row[2], row[2]+96)
        if len(row[4]) > 1008:
            # If the length of row[4] is greater than 25200 characters
            first_part = row[4][:498]  # First 12450 characters
            last_part = row[4][-498:]  # Last 12450 characters
            row4_modified = first_part + 'N' * 12 + last_part  # Insert 50 'N's in between
        else:
            # If the length of row[4] is less than or equal to 25200 characters
            row4_modified = row[4]
        seq = upstr + row4_modified + downstr
        if "N" in seq:
            seq = transform_genomic_sequence(seq)
        print(row[0], row[1], row[2], row[3], seq.upper(), sep="\t")
'''        
for i, row in input.iterrows():
    upstr = hg38fasta.fetch(row[0], row[1]-2400, row[1])
    downstr = hg38fasta.fetch(row[0], row[2], row[2]+2400)
    seq = upstr + downstr
    print(row[0], row[1], row[2], seq.upper(), sep="\t")
'''