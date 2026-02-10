#!/usr/bin/env python3
"""
Add PHRED-scaled columns to a scored CADD-SV output file (in-place).

For each score column (CADD-SV_score, CADD-SV-SR_score, CADD-SV-seqonly_score),
inserts a corresponding *_PHRED column immediately before it, using per-SV-type
lookup tables from the models/ directory.

Usage: add_phred.py <scored_file>
"""

import sys
import numpy as np
import pandas as pd
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent
MODEL_DIR = SCRIPT_DIR.parent / "models"

# Score columns eligible for PHRED scaling (as named in the PHRED tables)
_PHRED_SCORE_COLS = {"CADD-SV_score", "CADD-SV-SR_score", "CADD-SV-seqonly_score"}

# Output column name -> PHRED table column name (when they differ)
_SCORE_COL_ALIASES = {
    "CADD-SV_seqonly_score": "CADD-SV-seqonly_score",
}


def add_phred(scored_path):
    df = pd.read_table(scored_path, low_memory=False)

    if "type" not in df.columns:
        sys.stderr.write("[add_phred] No 'type' column found, skipping PHRED scaling.\n")
        return

    # Discover which score columns are present and eligible for PHRED
    score_cols = []
    for col in df.columns:
        phred_table_col = _SCORE_COL_ALIASES.get(col, col)
        if phred_table_col in _PHRED_SCORE_COLS:
            score_cols.append((col, phred_table_col))

    if not score_cols:
        sys.stderr.write("[add_phred] No score columns found, skipping PHRED scaling.\n")
        return

    # Load PHRED tables for each SV type present in the data, pre-sorted
    # by each score column so we don't re-sort inside the inner loop.
    phred_tables = {}
    for sv_type in df["type"].unique():
        phred_path = MODEL_DIR / f"{sv_type}_PHRED.tsv"
        if phred_path.exists():
            ptable = pd.read_table(phred_path)
            sorted_cache = {}
            for _, phred_table_col in score_cols:
                if phred_table_col in ptable.columns:
                    sorted_pt = ptable.sort_values(phred_table_col).reset_index(drop=True)
                    sorted_cache[phred_table_col] = (
                        sorted_pt[phred_table_col].values,
                        sorted_pt["PHRED"].values,
                    )
            if sorted_cache:
                phred_tables[sv_type] = sorted_cache

    if not phred_tables:
        sys.stderr.write("[add_phred] No PHRED tables found, skipping.\n")
        return

    # Pre-compute type masks once
    type_masks = {sv_type: df["type"] == sv_type for sv_type in phred_tables}

    # Build PHRED columns and insert right before each score column.
    # Process in reverse column-position order so inserts don't shift indices.
    inserts = []
    for output_col, phred_table_col in score_cols:
        col_idx = df.columns.get_loc(output_col)
        phred_col_name = output_col.replace("_score", "_PHRED")

        phred_values = pd.Series(np.nan, index=df.index)

        for sv_type, sorted_cache in phred_tables.items():
            if phred_table_col not in sorted_cache:
                continue

            mask = type_masks[sv_type]
            if not mask.any():
                continue

            sorted_scores, sorted_phred = sorted_cache[phred_table_col]
            raw_scores = df.loc[mask, output_col].values

            nan_mask = np.isnan(raw_scores.astype(float))
            safe_raw = np.where(nan_mask, 0.0, raw_scores.astype(float))

            indices = np.searchsorted(sorted_scores, safe_raw, side="right") - 1
            indices = np.clip(indices, 0, len(sorted_scores) - 1)

            result = sorted_phred[indices].astype(float)
            result[nan_mask] = np.nan
            result[~nan_mask & (safe_raw < sorted_scores[0])] = 0.0
            result[~nan_mask & (safe_raw > sorted_scores[-1])] = 50.0

            phred_values.loc[mask] = result

        inserts.append((col_idx, phred_col_name, phred_values))

    for col_idx, name, values in sorted(inserts, key=lambda x: x[0], reverse=True):
        df.insert(col_idx, name, values)

    df.to_csv(scored_path, sep="\t", index=False)
    sys.stderr.write(f"[add_phred] PHRED columns added to {scored_path}\n")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.stderr.write("Usage: add_phred.py <scored_file>\n")
        sys.exit(2)
    add_phred(sys.argv[1])
