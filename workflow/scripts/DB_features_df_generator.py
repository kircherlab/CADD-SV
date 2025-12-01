#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import h5py
import numpy as np

# Usage (REQUIRED): DB_features_df_generator.py coords.bed sample.h5 coords_ref.bed ref.h5
# Output columns:
#   chr, start, end, type, then per feature (14):
#     <feat>_absdelta_avg96, <feat>_absdelta_min96, <feat>_absdelta_max96
#
# Notes:
# - Coordinates file's sequence is already transformed by your generator,
#   so we use len(seq) directly.
# - Combined flank = union of first 96 and last 96 positions (no double-count).
# - We align sample/ref to a common length before indexing.

# -------- arg parsing --------
if len(sys.argv) != 5:
    sys.stderr.write(
        "ERROR: This script now requires reference inputs.\n"
        "Usage: DB_features_df_generator.py coords.bed sample.h5 coords_ref.bed ref.h5\n"
    )
    sys.exit(2)

coords_path = sys.argv[1]
sample_h5_path = sys.argv[2]
coords_ref_path = sys.argv[3]
ref_h5_path = sys.argv[4]

# -------- config --------
feature_names = [
    "protein_coding", "lncrna", "exon", "intron",
    "splice_donor", "splice_acceptor", "utr5", "utr3",
    "ctcf", "polya", "enh_tspec", "enh_tinv", "prom_tspec", "prom_tinv",
]
N_FEATS = len(feature_names)
FLANK = 96  # bp

def safe_mean(x: np.ndarray) -> float:
    return float(np.mean(x)) if x.size else 0.0

def safe_min(x: np.ndarray) -> float:
    return float(np.min(x)) if x.size else 0.0

def safe_max(x: np.ndarray) -> float:
    return float(np.max(x)) if x.size else 0.0

# -------- HDF5 helpers (layout-agnostic) --------
def _all_datasets(h5obj):
    if isinstance(h5obj, h5py.Dataset):
        yield ("", h5obj)
        return
    found = []
    def _vis(name, obj):
        if isinstance(obj, h5py.Dataset):
            found.append((name, obj))
    h5obj.visititems(_vis)
    for name, ds in found:
        yield (name, ds)

def _to_L14(arr: np.ndarray):
    a = np.asarray(arr)
    if a.ndim == 2 and a.shape[1] >= 14:
        return a[:, :14]
    while a.ndim > 2 and a.shape[0] == 1:
        a = a[0]
    while a.ndim > 2 and a.shape[-1] == 1:
        a = a[..., 0]
    if a.ndim == 2 and a.shape[1] >= 14:
        return a[:, :14]
    if a.ndim == 3:
        if a.shape[0] == 1 and a.shape[2] >= 14:  # (1, L, 14)
            return a[0, :, :14]
        if a.shape[2] == 1 and a.shape[1] >= 14:  # (L, 14, 1)
            return a[:, :14, 0]
    return None

def find_feature_matrix(h5group):
    candidates = []
    for path, ds in _all_datasets(h5group):
        try:
            shape = ds.shape
        except Exception:
            continue
        if len(shape) >= 2 and (shape[-1] >= 14 or (len(shape) >= 2 and shape[1] >= 14)):
            candidates.append((path, ds))
    for path, ds in sorted(candidates, key=lambda x: x[0]):
        try:
            arr = ds[...]
            mat = _to_L14(arr)
            if mat is not None:
                return mat
        except Exception as e:
            sys.stderr.write("[DB_features] Warning: failed reading dataset '{}': {}\n".format(path, e))
            continue
    return None

# -------- flank index helper (union without double counting) --------
def flank_indices(L: int, k: int = FLANK):
    if L <= 0:
        return np.array([], dtype=int)
    up_end = min(k, L)
    dn_start = max(0, L - k)
    idx = np.r_[np.arange(0, up_end), np.arange(dn_start, L)]
    return np.unique(idx) if idx.size else idx

# -------- headers --------
header = ["chr", "start", "end", "type"]
for feat in feature_names:
    header += [f"{feat}_absdelta_avg96", f"{feat}_absdelta_min96", f"{feat}_absdelta_max96"]
print("\t".join(header))

# -------- main --------
with h5py.File(sample_h5_path, "r") as f_sample, h5py.File(ref_h5_path, "r") as f_ref, \
     open(coords_path, "r") as fh_coords, open(coords_ref_path, "r") as fh_coords_ref:

    for line_idx, line in enumerate(fh_coords):
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            sys.stderr.write("[DB_features] Skipping line {}: expected >=5 columns, got {}\n".format(line_idx, len(parts)))
            continue

        chrom, start, end, type_, seq_s = parts[0], parts[1], parts[2], parts[3], parts[4]
        seq_len_s = len(seq_s)
        group_name = "group_{}".format(line_idx)

        # read matching ref line (keep alignment)
        ref_line = fh_coords_ref.readline()
        if not ref_line:
            sys.stderr.write("[DB_features] coords_ref ended early at line {}; skipping\n".format(line_idx))
            continue
        ref_parts = ref_line.rstrip("\n").split("\t")
        if len(ref_parts) < 5:
            sys.stderr.write("[DB_features] coords_ref line {} malformed; skipping\n".format(line_idx))
            continue
        seq_len_r = len(ref_parts[4])

        # fetch matrices
        if group_name not in f_sample or group_name not in f_ref:
            sys.stderr.write("[DB_features] Missing '{}' in one of the H5 files at line {}; skipping\n".format(group_name, line_idx))
            continue

        mat_s = find_feature_matrix(f_sample[group_name])
        mat_r = find_feature_matrix(f_ref[group_name])
        if mat_s is None or mat_r is None:
            sys.stderr.write("[DB_features] No usable (L,14) dataset in one/both groups at line {}; skipping\n".format(line_idx))
            continue

        # align to common length
        L_s = min(seq_len_s, mat_s.shape[0])
        L_r = min(seq_len_r, mat_r.shape[0])
        L = min(L_s, L_r)
        if L <= 0:
            sys.stderr.write("[DB_features] Zero common length at line {}; skipping\n".format(line_idx))
            continue

        mat_s = mat_s[:L, :]
        mat_r = mat_r[:L, :]

        # combined flank indices on common length
        idx = flank_indices(L, FLANK)
        if idx.size == 0:
            sys.stderr.write("[DB_features] Empty flank indices at line {}; skipping\n".format(line_idx))
            continue

        row_out = [chrom, start, end, type_]

        # per-feature abs(delta) stats
        for i in range(len(feature_names)):
            vals_s = mat_s[idx, i]
            vals_r = mat_r[idx, i]
            abs_delta = np.abs(vals_s - vals_r)
            row_out.append(str(safe_mean(abs_delta)))
            row_out.append(str(safe_min(abs_delta)))
            row_out.append(str(safe_max(abs_delta)))

        print("\t".join(row_out))
