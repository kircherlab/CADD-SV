#!/usr/bin/env python3

import sys
import pandas as pd
import numpy as np

# Reading command line arguments
args = sys.argv[1:]

# Load the input files
gn2 = pd.read_csv(args[0], sep="\t", header=None)
pli = pd.read_csv(args[1], sep=",", header=0)

gn = gn2.iloc[:, 3]  # Selecting the 4th column
gn = gn.astype(str)  # Convert to string
gene = pli.iloc[:, 0].values  # Convert the first column to a vector
pli_scores = pli.iloc[:, 1].values  # Extract PLI scores from the second column

# Gene to PLI score assignment
score = [None] * len(gn)

for i in np.where(gn != '0')[0]:  # Non-zero entries in the `gn` column
    a = gn[i].split(",")  # Split by comma
    score[i] = []
    for j in range(len(a)):
        gene_idx = np.where(gene == a[j])[0]
        if gene_idx.size > 0:
            score[i].append(pli_scores[gene_idx[0]])

# Initialize result vector
res = np.zeros(len(gn))

# Filter out non-empty scores and find the maximum score
if np.any([s for s in score if s]):
    for i in np.where([len(s) > 0 for s in score])[0]:
        if len(score[i]) > 0:
            res[i] = max(score[i])

# Identify empty entries (non-genic / no PLI)
a1 = np.where(res > 0)[0]
a2 = np.arange(len(gn))
a3 = np.setdiff1d(a2, a1)  # Vector of empty entries
res[a3] = 0  # Assign 0 to non-genic/no PLI entries

# Combine with gn2 excluding the 4th column and adding results
res2 = pd.concat([gn2.drop(columns=3), pd.Series(res)], axis=1)

# Save the result to the output file
res2.to_csv(args[2], sep="\t", index=False, header=False, quoting=3)
