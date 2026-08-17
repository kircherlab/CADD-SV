import pandas as pd
import numpy as np
import os
import sys

from row_order import read_inverse_permutation

# Function definition for scoring


def cadd_sv_read(file):
    y = pd.read_table(file, header=0, low_memory=False)
    y.iloc[:, 0] = y.iloc[:, 0].replace(to_replace='[chr\n]', value='', regex=True)
    y.iloc[:, 0] = y.iloc[:, 0].replace(to_replace='X', value='23', regex=True)
    y.iloc[:, 0] = y.iloc[:, 0].replace(to_replace='Y', value='24', regex=True)
    y = y.replace('.', 0)
    
    for i in range(1, y.shape[1]):
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
            "exon", "transcript", "gene", "3utr", "5utr", "cds", "nr_uc_bases",
            "RegSeq0_sum", "RegSeq1_sum", "RegSeq2_sum", "RegSeq3_sum",
            "RegSeq4_sum", "RegSeq5_sum", "RegSeq6_sum", "RegSeq7_sum",
            "RouletteAR_sum", "TADboundary_count", "boundary_score_sum", "screen_dELS",
            "screen_pELS", "screen_CA", "screen_CA-CTCF", "screen_CA-H3K4me3", "screen_CA-TF",
            "screen_TF", "screen_PLS", "screen_ELS_total", "screen_CA_total"]

    for k in tolog:  # logs of distances and sums and counts
        #i = header.index(k)
        y[k] = np.round(np.log10(np.abs(y[k]) + 1), 4)
    
    return y

def cadd_sv(matrix, up, down, genome, up_order, down_order):
    k = []
    y = []

    genome = pd.read_table(genome)  # added for ranges towards the end of chromosome
    genome.iloc[:, 0] = genome.iloc[:, 0].replace(to_replace='[chr\n]', value='', regex=True)  # added for ranges towards the end of chromosome
    
    k.append(cadd_sv_read(matrix))
    up_data = cadd_sv_read(up)
    down_data = cadd_sv_read(down)

    if len(up_data) != len(k[0]) or len(down_data) != len(k[0]):
        raise ValueError(
            "Whole-variant, upstream, and downstream matrices have different row counts."
        )

    up_inverse = read_inverse_permutation(up_order, len(up_data))
    down_inverse = read_inverse_permutation(down_order, len(down_data))
    up_data.iloc[:, 3:] = up_data.iloc[up_inverse, 3:].to_numpy()
    down_data.iloc[:, 3:] = down_data.iloc[down_inverse, 3:].to_numpy()

    k.append(up_data)
    k.append(down_data)

    k.append(k[1] + k[2])
    # Apply minimum for distance columns (values should not be summed across flanks)
    dist_cols = [col for col in k[3].columns if col.endswith('_dist') or col.endswith('_distance')]
    for col in dist_cols:
        k[3][col] = np.minimum(k[2][col], k[3][col])
    
    y.append(k[0])
    y.append(k[3])
 
    newy=k[0].join(k[3].iloc[:,3:], lsuffix="", rsuffix=("_flank"))
    newy.iloc[:, 0] = newy.iloc[:, 0].replace(to_replace='23', value='X', regex=True)
    newy.iloc[:, 0] = newy.iloc[:, 0].replace(to_replace='24', value='Y', regex=True)
    y.append(newy)    
    #y.append(cadd_sv_read(x, ext, z="flank"))
    
    return y
    

CB = cadd_sv(
    matrix=sys.argv[1],
    up=sys.argv[2],
    down=sys.argv[3],
    genome=sys.argv[4],
    up_order=sys.argv[5],
    down_order=sys.argv[6],
)

CB[2].to_csv(sys.argv[7], sep="\t", index=False)
