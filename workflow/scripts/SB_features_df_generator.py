import os
import sys
import h5py
import pandas as pd
import numpy as np
import scipy.signal as signal

# Path to the HDF5 file
h5_file_path = sys.argv[2]
tsv_file_path = sys.argv[1]
columns = [
    # Protein coding
    "chr", "start", "end", "type",
    "protein_coding_avg_whole", "protein_coding_avg_span",
    "protein_coding_min_whole", "protein_coding_min_span",
    "protein_coding_max_whole", "protein_coding_max_span",
    "protein_coding_count_whole", "protein_coding_count_span",

    # lncRNA
    "lncrna_avg_whole", "lncrna_avg_span",
    "lncrna_min_whole", "lncrna_min_span",
    "lncrna_max_whole", "lncrna_max_span",
    "lncrna_count_whole", "lncrna_count_span",

    # Exon
    "exon_avg_whole", "exon_avg_span",
    "exon_min_whole", "exon_min_span",
    "exon_max_whole", "exon_max_span",
    "exon_count_whole", "exon_count_span",

    # Intron
    "intron_avg_whole", "intron_avg_span",
    "intron_min_whole", "intron_min_span",
    "intron_max_whole", "intron_max_span",
    "intron_count_whole", "intron_count_span",

    # Splice donor
    "splice_donor_avg_whole", "splice_donor_avg_span",
    "splice_donor_min_whole", "splice_donor_min_span",
    "splice_donor_max_whole", "splice_donor_max_span",
    "splice_donor_count_whole", "splice_donor_count_span",

    # Splice acceptor
    "splice_acceptor_avg_whole", "splice_acceptor_avg_span",
    "splice_acceptor_min_whole", "splice_acceptor_min_span",
    "splice_acceptor_max_whole", "splice_acceptor_max_span",
    "splice_acceptor_count_whole", "splice_acceptor_count_span",

    # UTR5
    "utr5_avg_whole", "utr5_avg_span",
    "utr5_min_whole", "utr5_min_span",
    "utr5_max_whole", "utr5_max_span",
    "utr5_count_whole", "utr5_count_span",

    # UTR3
    "utr3_avg_whole", "utr3_avg_span",
    "utr3_min_whole", "utr3_min_span",
    "utr3_max_whole", "utr3_max_span",
    "utr3_count_whole", "utr3_count_span",

    # CTCF-bound
    "ctcf_avg_whole", "ctcf_avg_span",
    "ctcf_min_whole", "ctcf_min_span",
    "ctcf_max_whole", "ctcf_max_span",
    "ctcf_count_whole", "ctcf_count_span",

    # PolyA signal
    "polya_avg_whole", "polya_avg_span",
    "polya_min_whole", "polya_min_span",
    "polya_max_whole", "polya_max_span",
    "polya_count_whole", "polya_count_span",

    # Enhancer Tissue-specific
    "enh_tspec_avg_whole", "enh_tspec_avg_span",
    "enh_tspec_min_whole", "enh_tspec_min_span",
    "enh_tspec_max_whole", "enh_tspec_max_span",
    "enh_tspec_count_whole", "enh_tspec_count_span",

    # Enhancer Tissue-invariant
    "enh_tinv_avg_whole", "enh_tinv_avg_span",
    "enh_tinv_min_whole", "enh_tinv_min_span",
    "enh_tinv_max_whole", "enh_tinv_max_span",
    "enh_tinv_count_whole", "enh_tinv_count_span",

    # Promoter Tissue-specific
    "prom_tspec_avg_whole", "prom_tspec_avg_span",
    "prom_tspec_min_whole", "prom_tspec_min_span",
    "prom_tspec_max_whole", "prom_tspec_max_span",
    "prom_tspec_count_whole", "prom_tspec_count_span",

    # Promoter Tissue-invariant
    "prom_tinv_avg_whole", "prom_tinv_avg_span",
    "prom_tinv_min_whole", "prom_tinv_min_span",
    "prom_tinv_max_whole", "prom_tinv_max_span",
    "prom_tinv_count_whole", "prom_tinv_count_span"
]


distance_values = {
    0: 50,  # protein_coding_gene
    1: 50,  # lncRNA
    2: 50,   # exon
    3: 50,   # intron
    4: 1,    # splice_donor
    5: 1,    # splice_acceptor
    6: 50,   # 5UTR
    7: 50,   # 3UTR
    8: 10,   # CTCF-bound
    9: 1,    # polyA_signal
    10: 50,  # enhancer_Tissue_specific
    11: 50,  # enhancer_Tissue_invariant
    12: 50,  # promoter_Tissue_specific
    13: 50   # promoter_Tissue_invariant
}

width_values = {
    0: 50,   # protein_coding_gene
    1: 50,   # lncRNA
    2: 10,   # exon
    3: 10,   # intron
    4: 1,    # splice_donor
    5: 1,    # splice_acceptor
    6: 10,   # 5UTR
    7: 10,   # 3UTR
    8: 10,   # CTCF-bound
    9: 1,    # polyA_signal
    10: 10,  # enhancer_Tissue_specific
    11: 10,  # enhancer_Tissue_invariant
    12: 10,  # promoter_Tissue_specific
    13: 10   # promoter_Tissue_invariant
}


print("\t".join(columns))
# Open the HDF5 file
with h5py.File(h5_file_path, 'r') as f:
    # Open the TSV file
    with open(tsv_file_path, 'r') as tsv_file:
        for line_idx, line in enumerate(tsv_file):

            group_name = f'group_{line_idx}'
            line_split = line.strip().split('\t')
            last_column = line_split[4]

            def adjusted_length(s):
                return sum(6 if char == 'N' else 1 for char in s)
    
            # Apply the adjusted length function to each string in the last column
            string_lengths = adjusted_length(last_column)
            group = f[group_name]
            
            row_features = []
            row_features.append(line_split[0])
            row_features.append(line_split[1])
            row_features.append(line_split[2])
            row_features.append(line_split[3])
            for i_feature in range(0,14):
                data = [group[f'array_{j}'][...] for j in range(len(group.keys()))]
                
                data2 = data[0][0][:, i_feature]
                whole_array = data2[:string_lengths]

                # Split into span (middle section)
                span = whole_array[2400:-2399]
                distance_value = distance_values[i_feature]
                width_value = width_values[i_feature]
                peaks_whole, properties_whole = signal.find_peaks(whole_array, height=0.5, distance=distance_value, width=width_value)
                peaks_span, properties_span = signal.find_peaks(span, height=0.5, distance=distance_value, width=width_value)
                row_features.append(str(np.average(whole_array)))
                row_features.append(str(np.average(span)))
                row_features.append(str(np.min(whole_array)))
                row_features.append(str(np.min(span)))
                row_features.append(str(np.max(whole_array)))
                row_features.append(str(np.max(span)))
                row_features.append(str(len(peaks_whole)))
                row_features.append(str(len(peaks_span)))

            print("\t".join(row_features))