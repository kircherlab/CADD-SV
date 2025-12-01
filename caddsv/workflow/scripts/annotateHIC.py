#!/usr/bin/env python3

import sys
import numpy as np
import pandas as pd

# Reading command line arguments
args = sys.argv[1:]

# Load the input file
sub = pd.read_csv(args[0], sep="\t", header=None)

l = len(sub)

# Initialize vectors
x1 = np.zeros(l)
x2 = np.zeros(l)
x3 = np.zeros(l)
x4 = np.zeros(l)

# Populate x1, x2, x3, x4 based on the calculations
for i in range(l):
    x1[i] = sub.iloc[i, 4] - sub.iloc[i, 1]
    x2[i] = sub.iloc[i, 5] - sub.iloc[i, 2]
    x3[i] = sub.iloc[i, 4] - sub.iloc[i, 2]
    x4[i] = sub.iloc[i, 5] - sub.iloc[i, 1]

# Combine the x1, x2, x3, x4 arrays into a matrix
xx = np.column_stack((x1, x2, x3, x4))

# Find positions based on the conditions >0 and <0
p1 = np.where(xx[:, 0] > 0)[0]
p2 = np.where(xx[:, 1] > 0)[0]
p3 = np.where(xx[:, 2] > 0)[0]
p4 = np.where(xx[:, 3] > 0)[0]

n1 = np.where(xx[:, 0] < 0)[0]
n2 = np.where(xx[:, 1] < 0)[0]
n3 = np.where(xx[:, 2] < 0)[0]
n4 = np.where(xx[:, 3] < 0)[0]

# Compute the various TAD categories using intersections
btad = np.intersect1d(np.intersect1d(np.intersect1d(p1, p2), p3), p4)  # Before a TAD
atad = np.intersect1d(np.intersect1d(np.intersect1d(n1, n2), n3), n4)  # After a TAD
stad = np.intersect1d(np.intersect1d(np.intersect1d(n1, p2), n3), p4)  # Spanning a TAD
itad = np.intersect1d(np.intersect1d(np.intersect1d(p1, n2), n3), p4)  # Intra a TAD
lbtad = np.intersect1d(np.intersect1d(np.intersect1d(p1, p2), n3), p4) # Left boundary a TAD
rbtad = np.intersect1d(np.intersect1d(np.intersect1d(n1, n2), n3), p4) # Right boundary a TAD

# Initialize vectors for in TAD, out TAD, boundary TAD, and no TAD
intad = np.zeros(l, dtype=int)
intad[stad] = 1
intad[itad] = 1

notad = np.zeros(l, dtype=int)
notad[btad] = 1
notad[atad] = 1

boundarytad = np.zeros(l, dtype=int)
boundarytad[lbtad] = 1
boundarytad[rbtad] = 1

# Calculate the minimum absolute distance to the boundary for each row
m = np.zeros(l)
for i in range(l):
    m[i] = np.min(np.abs(xx[i, :]))

# Combine the results into a dataframe (chr, start, end, in_tad, out_tad, boundary, distance)
bed = pd.DataFrame({
    'chr': sub.iloc[:, 0].astype(str),
    'start': sub.iloc[:, 1],
    'end': sub.iloc[:, 2],
    'in_tad': intad,
    'out_tad': notad,
    'boundary': boundarytad,
    'distance': m
})

# Save the result to the output file
bed.to_csv(args[1], sep="\t", header=False, index=False, quoting=3)
