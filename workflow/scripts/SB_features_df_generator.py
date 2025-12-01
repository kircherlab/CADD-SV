#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import h5py
import numpy as np

# Usage: SB_features_df_generator.py coords.bed sample.h5
# Output columns: SB-style avg/min/max for whole & span for 14 features.

if len(sys.argv) != 3:
    sys.stderr.write(
        "Usage: SB_features_df_generator.py coords.bed sample.h5\n"
    )
    sys.exit(2)

tsv_file_path = sys.argv[1]
h5_file_path = sys.argv[2]

# --- config ---
feature_names = [
    "protein_coding", "lncrna", "exon", "intron",
    "splice_donor", "splice_acceptor", "utr5", "utr3",
    "ctcf", "polya", "enh_tspec", "enh_tinv", "prom_tspec", "prom_tinv",
]
FLANK = 96  # bp for upstream; span removes upstream_length and 95 from the end to match original SB code

# SB column order
columns = [
    "chr", "start", "end", "type",
]
for feat in feature_names:
    columns += [
        f"{feat}_avg_whole", f"{feat}_avg_span",
        f"{feat}_min_whole", f"{feat}_min_span",
        f"{feat}_max_whole", f"{feat}_max_span",
    ]

def safe_mean(x: np.ndarray) -> float:
    return float(np.mean(x)) if x.size else 0.0

def safe_min(x: np.ndarray) -> float:
    return float(np.min(x)) if x.size else 0.0

def safe_max(x: np.ndarray) -> float:
    return float(np.max(x)) if x.size else 0.0

# ---- HDF5 helpers (DB-style, layout-agnostic) ----
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
    # Try common shapes first
    if a.ndim == 2 and a.shape[1] >= 14:
        return a[:, :14]
    # Squeeze leading singleton dims
    while a.ndim > 2 and a.shape[0] == 1:
        a = a[0]
    while a.ndim > 2 and a.shape[-1] == 1:
        a = a[..., 0]
    if a.ndim == 2 and a.shape[1] >= 14:
        return a[:, :14]
    # Handle (1, L, C) or (L, C, 1)
    if a.ndim == 3:
        if a.shape[0] == 1 and a.shape[2] >= 14:      # (1, L, C)
            return a[0, :, :14]
        if a.shape[2] == 1 and a.shape[1] >= 14:      # (L, C, 1)
            return a[:, :14, 0]
    return None

def find_feature_matrix(h5group):
    """
    Find a dataset inside `h5group` that can be coerced to (L,14).
    Returns a NumPy array of shape (L,14), or None.
    """
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
            sys.stderr.write(f"[SB_features] Warning: failed reading dataset '{path}': {e}\n")
            continue
    return None

def adjusted_length(seq: str) -> int:
    """
    SB behavior: count 'N' as length 6, others as 1.
    """
    return sum(6 if c == 'N' else 1 for c in seq)

# ---- run ----
print("\t".join(columns))

with h5py.File(h5_file_path, "r") as f, open(tsv_file_path, "r") as tsv_file:
    for line_idx, line in enumerate(tsv_file):
        line = line.rstrip("\n")
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 5:
            sys.stderr.write(f"[SB_features] Skipping line {line_idx}: expected >=5 columns, got {len(parts)}\n")
            continue

        chrom, start_s, end_s, type_s, seq_s = parts[0], parts[1], parts[2], parts[3], parts[4]

        # upstream length like original SB: min(96, start-1) but not below 0
        try:
            start_int = int(start_s)
        except ValueError:
            start_int = 0
        upstream_length = FLANK if start_int > FLANK else max(0, start_int - 1)

        # sequence-derived length with SB's adjusted length rule
        L_seq = adjusted_length(seq_s)

        group_name = f"group_{line_idx}"
        if group_name not in f:
            sys.stderr.write(f"[SB_features] Missing '{group_name}' in H5 at line {line_idx}; skipping\n")
            continue

        mat = find_feature_matrix(f[group_name])
        if mat is None:
            sys.stderr.write(f"[SB_features] No usable (L,14) dataset in '{group_name}' at line {line_idx}; skipping\n")
            continue

        # Align to common length with sequence
        L = min(L_seq, mat.shape[0])
        if L <= 0:
            sys.stderr.write(f"[SB_features] Zero length at line {line_idx}; skipping\n")
            continue

        mat = mat[:L, :]

        # Whole region is 0..L
        # Span region follows your SB logic: remove 'upstream_length' from start and 95 from end
        span_start = min(upstream_length, L)  # cap within [0, L]
        span_end = max(0, L - 95)            # as in your original code
        # If span_end <= span_start, span will be empty; safe_* will return 0s.

        row = [chrom, start_s, end_s, type_s]

        for i in range(len(feature_names)):
            whole = mat[:, i]
            span = whole[span_start:span_end]

            row.append(str(safe_mean(whole)))
            row.append(str(safe_mean(span)))

            row.append(str(safe_min(whole)))
            row.append(str(safe_min(span)))

            row.append(str(safe_max(whole)))
            row.append(str(safe_max(span)))

        print("\t".join(row))
