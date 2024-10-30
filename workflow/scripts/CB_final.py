import pandas as pd
import numpy as np
import os
import sys

# Function definition for scoring


def cadd_sv_read(x, ext, z):
    y = pd.read_table(f"{x}/{x}{ext}_matrix.bed{z}", header=0, low_memory=False)
    y.iloc[:, 0] = y.iloc[:, 0].replace(to_replace='[chr\n]', value='', regex=True)
    y.iloc[:, 0] = y.iloc[:, 0].replace(to_replace='X', value='23', regex=True)
    y.iloc[:, 0] = y.iloc[:, 0].replace(to_replace='Y', value='24', regex=True)
    y = y.replace('.', 0)
    
    for i in range(1, 131):
        y.iloc[:, i] = pd.to_numeric(y.iloc[:, i], errors='coerce')
    
    y = y.fillna(0)
    y = y.apply(pd.to_numeric)
    y = y.fillna(0)
    header = list(y.columns)
    y = pd.DataFrame(y)
    
    
#   ylen = y.iloc[:, 2] - y.iloc[:, 1]
    y.loc[y['DNase-seq_max'] > 3000, 'DNase-seq_max'] = 3000  # DNAse outliers
    y.loc[y['DDD_HaploInsuf'] > 1, 'DDD_HaploInsuf'] = 1  # DDD outliers
    
    tolog = ["EP_distance", "A549_nested_dist", "A549_tad_dist",
            "caki2_nested_dist", "caki2_tad_dist", "escTAD_distance",
            "microsyn_distance", "exon_dist", "gene_dist",
            "start_codon_dist", "CADD_sum", "CADD_count", "gerp_count", "PhastCons100_sum",
            "PhastCons30_sum", "PhastCons20_sum", "DI_min",
            "DI_max", "DNase-seq_sum", "H2AFZ_sum",
            "H3K27ac_sum", "H3K27me3_sum", "H3k36me3_sum",
            "H3K4me1_sum", "H3K4me2_sum", "H3K4me3_sum",
            "H3K79me2_sum", "H3K9ac_sum", "H3K9me3_sum",
            "H4k20me1_sum", "totalRNA-seq_sum", "LINSIGHT",
            "exon", "transcript", "gene", "3utr", "5utr", "cds", "nr_uc_bases"]

    for k in tolog:  # logs of distances and sums and counts
        #i = header.index(k)
        y[k] = np.round(np.log10(np.abs(y[k]) + 1), 4)
    
    return y

def cadd_sv(x, ext, genome):
    k = []
    y = []

    genome = pd.read_table(genome)  # added for ranges towards the end of chromosome
    genome.iloc[:, 0] = genome.iloc[:, 0].replace(to_replace='[chr\n]', value='', regex=True)  # added for ranges towards the end of chromosome
    
    k.append(cadd_sv_read(x, ext, z=""))
    k.append(cadd_sv_read(x, ext, z="up"))
    k.append(cadd_sv_read(x, ext, z="down"))

    k.append(k[1] + k[2])
    k[3].iloc[:, 121] = np.minimum(k[2].iloc[:, 121], k[3].iloc[:, 121])
    k[3].iloc[:, 122] = np.minimum(k[2].iloc[:, 122], k[3].iloc[:, 122])
    k[3].iloc[:, 123] = np.minimum(k[2].iloc[:, 123], k[3].iloc[:, 123])
    
    y.append(k[0])
    y.append(k[3])
 
    newy=k[0].join(k[3].iloc[:,3:131], lsuffix="", rsuffix=("_flank"))
    newy.iloc[:, 0] = newy.iloc[:, 0].replace(to_replace='23', value='X', regex=True)
    newy.iloc[:, 0] = newy.iloc[:, 0].replace(to_replace='24', value='Y', regex=True)
    y.append(newy)    
    #y.append(cadd_sv_read(x, ext, z="flank"))
    
    return y
    

CB = cadd_sv(x=sys.argv[1], ext="bed", genome="ucsc/hg38.fa.genome")

CB[2].to_csv(sys.argv[2], sep="\t", index=False)
