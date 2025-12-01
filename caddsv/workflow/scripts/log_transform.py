import pandas as pd
import numpy as np
import sys

input_file = sys.argv[1]
output_file = sys.argv[2]

df = pd.read_csv(input_file, sep='\t')

# Identify columns with "_count_" in header
tolog = [col for col in df.columns if "_count_" in col]

# Cap values at 50, then apply log transformation
for k in tolog:
    df[k] = np.clip(df[k], None, 50)  # Cap at 50
    df[k] = np.round(np.log10(np.abs(df[k]) + 1), 4)  # Log-transform

df.to_csv(output_file, sep='\t', index=False)
