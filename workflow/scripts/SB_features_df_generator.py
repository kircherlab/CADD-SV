import os
import sys
import h5py
import pandas as pd
import numpy as np
import scipy.signal as signal

def load_sb_file(file_name):
    """Loads the SB_probabilities.h5 file using h5py."""
    loaded_probabilities = []
    with h5py.File(file_name, 'r') as f:
        # Iterate over each group
        for i in range(len(f.keys())):
            grp = f[f'group_{i}']
            sublist = [grp[f'array_{j}'][...] for j in range(len(grp.keys()))]
            loaded_probabilities.append(sublist)
    return sublist

def load_non_sb_file(file_name):
    """Loads the non-SB file using pandas with no header and tab-separated."""
    df = pd.read_csv(file_name, sep='\t', header=None)
    return df

def count_string_length(non_sb_data):
    """Counts the length of the string in the last column for each row."""
    # Extract the last column (assuming it's a string column)
    last_column = non_sb_data.iloc[:, -1]
    
    # Count the length of the string for each row
    string_lengths = last_column.apply(len)
    
    return string_lengths

def main(folder_path):
    # Get all files ending with SB_probabilities.h5
    files = [f for f in os.listdir(folder_path) if f.endswith('SB_probabilities.h5')]
    #larger_df = None
    columns = [
        "chr", "start", "end", "protein_coding_avg_whole", "protein_coding_min_whole", 
        "protein_coding_max_whole", "protein_coding_count_whole", "lncrna_avg_whole", 
        "lncrna_min_whole", "lncrna_max_whole", "lncrna_count_whole", "utr5_avg_whole", 
        "utr5_min_whole", "utr5_max_whole", "utr5_count_whole", "utr3_avg_whole", 
        "utr3_min_whole", "utr3_max_whole", "utr3_count_whole", "exon_avg_whole", 
        "exon_min_whole", "exon_max_whole", "exon_count_whole", "intron_avg_whole", 
        "intron_min_whole", "intron_max_whole", "intron_count_whole", 
        "splice_donor_avg_whole", "splice_donor_min_whole", "splice_donor_max_whole", 
        "splice_donor_count_whole", "splice_acceptor_avg_whole", "splice_acceptor_min_whole", 
        "splice_acceptor_max_whole", "splice_acceptor_count_whole", "prom_tspec_avg_whole", 
        "prom_tspec_min_whole", "prom_tspec_max_whole", "prom_tspec_count_whole", 
        "prom_tinv_avg_whole", "prom_tinv_min_whole", "prom_tinv_max_whole", 
        "prom_tinv_count_whole", "enh_tspec_avg_whole", "enh_tspec_min_whole", 
        "enh_tspec_max_whole", "enh_tspec_count_whole", "enh_tinv_avg_whole", 
        "enh_tinv_min_whole", "enh_tinv_max_whole", "enh_tinv_count_whole", 
        "polya_avg_whole", "polya_min_whole", "polya_max_whole", "polya_count_whole", 
        "ctcf_avg_whole", "ctcf_min_whole", "ctcf_max_whole", "ctcf_count_whole", 
        "protein_coding_avg_span", "protein_coding_min_span", "protein_coding_max_span", 
        "protein_coding_count_span", "lncrna_avg_span", "lncrna_min_span", "lncrna_max_span", 
        "lncrna_count_span", "utr5_avg_span", "utr5_min_span", "utr5_max_span", 
        "utr5_count_span", "utr3_avg_span", "utr3_min_span", "utr3_max_span", 
        "utr3_count_span", "exon_avg_span", "exon_min_span", "exon_max_span", 
        "exon_count_span", "intron_avg_span", "intron_min_span", "intron_max_span", 
        "intron_count_span", "splice_donor_avg_span", "splice_donor_min_span", 
        "splice_donor_max_span", "splice_donor_count_span", "splice_acceptor_avg_span", 
        "splice_acceptor_min_span", "splice_acceptor_max_span", "splice_acceptor_count_span", 
        "prom_tspec_avg_span", "prom_tspec_min_span", "prom_tspec_max_span", 
        "prom_tspec_count_span", "prom_tinv_avg_span", "prom_tinv_min_span", 
        "prom_tinv_max_span", "prom_tinv_count_span", "enh_tspec_avg_span", 
        "enh_tspec_min_span", "enh_tspec_max_span", "enh_tspec_count_span", 
        "enh_tinv_avg_span", "enh_tinv_min_span", "enh_tinv_max_span", 
        "enh_tinv_count_span", "polya_avg_span", "polya_min_span", "polya_max_span", 
        "polya_count_span", "ctcf_avg_span", "ctcf_min_span", "ctcf_max_span", 
        "ctcf_count_span"
    ]

    print("\t".join(columns))
    # Iterate through all SB files
    for sb_file in files:
        # Extract the prefix (i.e., the "part_xxx" part of the filename)
        prefix = sb_file.replace('SB_probabilities.h5', '')
        
        # Look for the corresponding non-SB file
        non_sb_file = os.path.join(folder_path, prefix)
        
        if os.path.exists(non_sb_file):
            # Load the SB file
            sb_file_path = os.path.join(folder_path, sb_file)
            sb_data = load_sb_file(sb_file_path)
            
            # Load the non-SB file
            non_sb_data = load_non_sb_file(non_sb_file)
            
            # Step 2: Count the length of the string in the last column for each row
            string_lengths = count_string_length(non_sb_data)
            
            # Initialize empty lists for each feature category
            protein_coding_whole = []
            protein_coding_span = []
            lncrna_whole = []
            lncrna_span = []
            utr5_whole = []
            utr5_span = []
            utr3_whole = []
            utr3_span = []
            exon_whole = []
            exon_span = []
            intron_whole = []
            intron_span = []
            splice_donor_whole = []
            splice_donor_span = []
            splice_acceptor_whole = []
            splice_acceptor_span = []
            prom_tspec_whole = []
            prom_tspec_span = []
            prom_tinv_whole = []
            prom_tinv_span = []
            enh_tspec_whole = []
            enh_tspec_span = []
            enh_tinv_whole = []
            enh_tinv_span = []
            polya_whole = []
            polya_span = []
            ctcf_whole = []
            ctcf_span = []

            
            
            
            for i_event in range(len(non_sb_data)):
                for i_feature in range(0,14):
                    probabilities_feature = sb_data[0][..., i_feature][i_event]
                    # Define the distance and width for each value of i_feature (0 to 13)
                    distance_values = {
                        0: 200, 1: 200, 2: 50, 3: 50, 4: 50, 
                        5: 50, 6: 1, 7: 1, 8: 50, 
                        9: 50, 10: 50, 11: 50, 12: 10, 13: 1
                    }

                    width_values = {
                        0: 50, 1: 50, 2: 10, 3: 10, 4: 10, 
                        5: 10, 6: 1, 7: 1, 8: 10, 
                        9: 10, 10: 10, 11: 10, 12: 10, 13: 1
                    }

                    # Subset the array based on string_lengths[i_event]
                    whole_array = probabilities_feature[:string_lengths[i_event]]

                    # Split into span (middle section)
                    span = whole_array[2400:-2399]

                    # Score the entire array (whole_array)
                    f_avg_whole = float(np.average(whole_array))
                    f_min_whole = float(np.min(whole_array))
                    f_max_whole = float(np.max(whole_array))

                    # Score span (the middle part)
                    f_avg_span = float(np.average(span))
                    f_min_span = float(np.min(span))
                    f_max_span = float(np.max(span))

                    # Get the specific distance and width values for the current i_feature
                    distance_value = distance_values[i_feature]
                    width_value = width_values[i_feature]

                    # Find peaks in the whole_array using the specific distance and width for this i_feature
                    peaks_whole, properties_whole = signal.find_peaks(whole_array, height=0.5, distance=distance_value, width=width_value)
                    f_count_whole = len(peaks_whole)  # Number of peaks in whole_array

                    # Find peaks in span using the specific distance and width for this i_feature
                    peaks_span, properties_span = signal.find_peaks(span, height=0.5, distance=distance_value, width=width_value)
                    f_count_span = len(peaks_span)  # Number of peaks in span
                    # List of features for whole and span
                    list_features_whole = [f_avg_whole, f_min_whole, f_max_whole, f_count_whole]
                    list_features_span = [f_avg_span, f_min_span, f_max_span, f_count_span]

                    # Assign the features based on i_feature
                    if i_feature == 0:
                        protein_coding_whole.append(list_features_whole)
                        protein_coding_span.append(list_features_span)
                    elif i_feature == 1:
                        lncrna_whole.append(list_features_whole)
                        lncrna_span.append(list_features_span)
                    elif i_feature == 2:
                        utr5_whole.append(list_features_whole)
                        utr5_span.append(list_features_span)
                    elif i_feature == 3:
                        utr3_whole.append(list_features_whole)
                        utr3_span.append(list_features_span)
                    elif i_feature == 4:
                        exon_whole.append(list_features_whole)
                        exon_span.append(list_features_span)
                    elif i_feature == 5:
                        intron_whole.append(list_features_whole)
                        intron_span.append(list_features_span)
                    elif i_feature == 6:
                        splice_donor_whole.append(list_features_whole)
                        splice_donor_span.append(list_features_span)
                    elif i_feature == 7:
                        splice_acceptor_whole.append(list_features_whole)
                        splice_acceptor_span.append(list_features_span)
                    elif i_feature == 8:
                        prom_tspec_whole.append(list_features_whole)
                        prom_tspec_span.append(list_features_span)
                    elif i_feature == 9:
                        prom_tinv_whole.append(list_features_whole)
                        prom_tinv_span.append(list_features_span)
                    elif i_feature == 10:
                        enh_tspec_whole.append(list_features_whole)
                        enh_tspec_span.append(list_features_span)
                    elif i_feature == 11:
                        enh_tinv_whole.append(list_features_whole)
                        enh_tinv_span.append(list_features_span)
                    elif i_feature == 12:
                        polya_whole.append(list_features_whole)
                        polya_span.append(list_features_span)
                    elif i_feature == 13:
                        ctcf_whole.append(list_features_whole)
                        ctcf_span.append(list_features_span)

            # Convert the lists to DataFrames
            protein_coding_whole = pd.DataFrame(protein_coding_whole)
            lncrna_whole = pd.DataFrame(lncrna_whole)
            utr5_whole = pd.DataFrame(utr5_whole)
            utr3_whole = pd.DataFrame(utr3_whole)
            exon_whole = pd.DataFrame(exon_whole)
            intron_whole = pd.DataFrame(intron_whole)
            splice_donor_whole = pd.DataFrame(splice_donor_whole)
            splice_acceptor_whole = pd.DataFrame(splice_acceptor_whole)
            prom_tspec_whole = pd.DataFrame(prom_tspec_whole)
            prom_tinv_whole = pd.DataFrame(prom_tinv_whole)
            enh_tspec_whole = pd.DataFrame(enh_tspec_whole)
            enh_tinv_whole = pd.DataFrame(enh_tinv_whole)
            polya_whole = pd.DataFrame(polya_whole)
            ctcf_whole = pd.DataFrame(ctcf_whole)

            # Do the same for span
            protein_coding_span = pd.DataFrame(protein_coding_span)
            lncrna_span = pd.DataFrame(lncrna_span)
            utr5_span = pd.DataFrame(utr5_span)
            utr3_span = pd.DataFrame(utr3_span)
            exon_span = pd.DataFrame(exon_span)
            intron_span = pd.DataFrame(intron_span)
            splice_donor_span = pd.DataFrame(splice_donor_span)
            splice_acceptor_span = pd.DataFrame(splice_acceptor_span)
            prom_tspec_span = pd.DataFrame(prom_tspec_span)
            prom_tinv_span = pd.DataFrame(prom_tinv_span)
            enh_tspec_span = pd.DataFrame(enh_tspec_span)
            enh_tinv_span = pd.DataFrame(enh_tinv_span)
            polya_span = pd.DataFrame(polya_span)
            ctcf_span = pd.DataFrame(ctcf_span)
            chr_column = non_sb_data.iloc[:, 0]  # First column: chr
            start_column = non_sb_data.iloc[:, 1]  # Second column: start
            end_column = non_sb_data.iloc[:, 2]  # Third column: end
            # New features dictionary combining whole and span features
            new_features = {
                "chr": chr_column,
                "start": start_column,
                "end": end_column,
                "protein_coding_avg_whole": protein_coding_whole[0],
                "protein_coding_min_whole": protein_coding_whole[1],
                "protein_coding_max_whole": protein_coding_whole[2],
                "protein_coding_count_whole": protein_coding_whole[3],

                "lncrna_avg_whole": lncrna_whole[0],
                "lncrna_min_whole": lncrna_whole[1],
                "lncrna_max_whole": lncrna_whole[2],
                "lncrna_count_whole": lncrna_whole[3],

                "utr5_avg_whole": utr5_whole[0],
                "utr5_min_whole": utr5_whole[1],
                "utr5_max_whole": utr5_whole[2],
                "utr5_count_whole": utr5_whole[3],

                "utr3_avg_whole": utr3_whole[0],
                "utr3_min_whole": utr3_whole[1],
                "utr3_max_whole": utr3_whole[2],
                "utr3_count_whole": utr3_whole[3],

                "exon_avg_whole": exon_whole[0],
                "exon_min_whole": exon_whole[1],
                "exon_max_whole": exon_whole[2],
                "exon_count_whole": exon_whole[3],

                "intron_avg_whole": intron_whole[0],
                "intron_min_whole": intron_whole[1],
                "intron_max_whole": intron_whole[2],
                "intron_count_whole": intron_whole[3],

                "splice_donor_avg_whole": splice_donor_whole[0],
                "splice_donor_min_whole": splice_donor_whole[1],
                "splice_donor_max_whole": splice_donor_whole[2],
                "splice_donor_count_whole": splice_donor_whole[3],

                "splice_acceptor_avg_whole": splice_acceptor_whole[0],
                "splice_acceptor_min_whole": splice_acceptor_whole[1],
                "splice_acceptor_max_whole": splice_acceptor_whole[2],
                "splice_acceptor_count_whole": splice_acceptor_whole[3],

                "prom_tspec_avg_whole": prom_tspec_whole[0],
                "prom_tspec_min_whole": prom_tspec_whole[1],
                "prom_tspec_max_whole": prom_tspec_whole[2],
                "prom_tspec_count_whole": prom_tspec_whole[3],

                "prom_tinv_avg_whole": prom_tinv_whole[0],
                "prom_tinv_min_whole": prom_tinv_whole[1],
                "prom_tinv_max_whole": prom_tinv_whole[2],
                "prom_tinv_count_whole": prom_tinv_whole[3],

                "enh_tspec_avg_whole": enh_tspec_whole[0],
                "enh_tspec_min_whole": enh_tspec_whole[1],
                "enh_tspec_max_whole": enh_tspec_whole[2],
                "enh_tspec_count_whole": enh_tspec_whole[3],

                "enh_tinv_avg_whole": enh_tinv_whole[0],
                "enh_tinv_min_whole": enh_tinv_whole[1],
                "enh_tinv_max_whole": enh_tinv_whole[2],
                "enh_tinv_count_whole": enh_tinv_whole[3],

                "polya_avg_whole": polya_whole[0],
                "polya_min_whole": polya_whole[1],
                "polya_max_whole": polya_whole[2],
                "polya_count_whole": polya_whole[3],

                "ctcf_avg_whole": ctcf_whole[0],
                "ctcf_min_whole": ctcf_whole[1],
                "ctcf_max_whole": ctcf_whole[2],
                "ctcf_count_whole": ctcf_whole[3],

                # Span features
                "protein_coding_avg_span": protein_coding_span[0],
                "protein_coding_min_span": protein_coding_span[1],
                "protein_coding_max_span": protein_coding_span[2],
                "protein_coding_count_span": protein_coding_span[3],

                "lncrna_avg_span": lncrna_span[0],
                "lncrna_min_span": lncrna_span[1],
                "lncrna_max_span": lncrna_span[2],
                "lncrna_count_span": lncrna_span[3],

                "utr5_avg_span": utr5_span[0],
                "utr5_min_span": utr5_span[1],
                "utr5_max_span": utr5_span[2],
                "utr5_count_span": utr5_span[3],

                "utr3_avg_span": utr3_span[0],
                "utr3_min_span": utr3_span[1],
                "utr3_max_span": utr3_span[2],
                "utr3_count_span": utr3_span[3],

                "exon_avg_span": exon_span[0],
                "exon_min_span": exon_span[1],
                "exon_max_span": exon_span[2],
                "exon_count_span": exon_span[3],

                "intron_avg_span": intron_span[0],
                "intron_min_span": intron_span[1],
                "intron_max_span": intron_span[2],
                "intron_count_span": intron_span[3],

                "splice_donor_avg_span": splice_donor_span[0],
                "splice_donor_min_span": splice_donor_span[1],
                "splice_donor_max_span": splice_donor_span[2],
                "splice_donor_count_span": splice_donor_span[3],

                "splice_acceptor_avg_span": splice_acceptor_span[0],
                "splice_acceptor_min_span": splice_acceptor_span[1],
                "splice_acceptor_max_span": splice_acceptor_span[2],
                "splice_acceptor_count_span": splice_acceptor_span[3],

                "prom_tspec_avg_span": prom_tspec_span[0],
                "prom_tspec_min_span": prom_tspec_span[1],
                "prom_tspec_max_span": prom_tspec_span[2],
                "prom_tspec_count_span": prom_tspec_span[3],

                "prom_tinv_avg_span": prom_tinv_span[0],
                "prom_tinv_min_span": prom_tinv_span[1],
                "prom_tinv_max_span": prom_tinv_span[2],
                "prom_tinv_count_span": prom_tinv_span[3],

                "enh_tspec_avg_span": enh_tspec_span[0],
                "enh_tspec_min_span": enh_tspec_span[1],
                "enh_tspec_max_span": enh_tspec_span[2],
                "enh_tspec_count_span": enh_tspec_span[3],

                "enh_tinv_avg_span": enh_tinv_span[0],
                "enh_tinv_min_span": enh_tinv_span[1],
                "enh_tinv_max_span": enh_tinv_span[2],
                "enh_tinv_count_span": enh_tinv_span[3],

                "polya_avg_span": polya_span[0],
                "polya_min_span": polya_span[1],
                "polya_max_span": polya_span[2],
                "polya_count_span": polya_span[3],

                "ctcf_avg_span": ctcf_span[0],
                "ctcf_min_span": ctcf_span[1],
                "ctcf_max_span": ctcf_span[2],
                "ctcf_count_span": ctcf_span[3]
            }

            new_df = pd.DataFrame(new_features)
            '''
            if larger_df is None:
                # First iteration: Initialize the larger DataFrame
                larger_df = new_df
            else:
                # Subsequent iterations: Append the new DataFrame to the larger DataFrame
                larger_df = pd.concat([larger_df, new_df], ignore_index=True)
            '''
            print(new_df.to_csv(sep='\t', index=False, header=False), end='')
        else:
            print(f"Warning: Corresponding non-SB file not found for {sb_file}")
    #larger_df.to_csv((folder_path + "SB_features_final.tsv"), sep="\t", index=False)
    
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print("Usage: python script.py <folder_path>")
        sys.exit(1)

    folder_path = sys.argv[1]
    
    if not os.path.isdir(folder_path):
        print(f"The provided path is not a valid directory: {folder_path}")
        sys.exit(1)
    
    main(folder_path)
