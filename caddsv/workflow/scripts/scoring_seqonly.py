#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scoring script for sequence-only mode.

Combines SB, SBref, and DB features and runs the seqonlymodel.joblib model.

Usage: scoring_seqonly.py SB.bed SBref.bed DB.bed output.tsv
"""

import sys
import numpy as np
import pandas as pd
import joblib
from pathlib import Path

# Paths
SCRIPT_DIR = Path(__file__).resolve().parent
WORKFLOW_DIR = SCRIPT_DIR.parent
MODEL_DIR = WORKFLOW_DIR / "models"

model_path = MODEL_DIR / "seqonlymodel.joblib"

if len(sys.argv) != 5:
    sys.stderr.write(
        "Usage: scoring_seqonly.py SB.bed SBref.bed DB.bed output.tsv\n"
    )
    sys.exit(2)

sb_file = sys.argv[1]
sbref_file = sys.argv[2]
db_file = sys.argv[3]
output_file = sys.argv[4]

# Load model
if not model_path.exists():
    sys.stderr.write(f"Error: Model not found at {model_path}\n")
    sys.exit(1)

model = joblib.load(model_path)

# Read feature files
sb_df = pd.read_table(sb_file, low_memory=False)
sbref_df = pd.read_table(sbref_file, low_memory=False)
db_df = pd.read_table(db_file, low_memory=False)

# Extract ID columns from SB (they should be the same across all files)
id_cols = ["chr", "start", "end", "type"]
ids = sb_df[id_cols].copy()

# Get feature columns (everything except the ID columns)
sb_features = sb_df.drop(columns=id_cols)
sbref_features = sbref_df.drop(columns=id_cols)
db_features = db_df.drop(columns=id_cols)

# Rename columns to avoid collisions
sb_features.columns = [f"SB_{c}" for c in sb_features.columns]
sbref_features.columns = [f"SBref_{c}" for c in sbref_features.columns]
db_features.columns = [f"DB_{c}" for c in db_features.columns]

# Combine all features
combined = pd.concat([sb_features, sbref_features, db_features], axis=1)

# Select features the model was trained on
try:
    X = combined[model.feature_names_in_]
except KeyError as e:
    # If some features are missing, report which ones
    missing = set(model.feature_names_in_) - set(combined.columns)
    extra = set(combined.columns) - set(model.feature_names_in_)
    sys.stderr.write(f"Feature mismatch!\n")
    sys.stderr.write(f"Missing features: {missing}\n")
    sys.stderr.write(f"Extra features: {extra}\n")
    sys.stderr.write(f"Available features: {list(combined.columns)}\n")
    sys.exit(1)

# Run predictions
predictions = model.predict_proba(X)[:, 1]

# Create output dataframe
output = ids.copy()
output["CADD-SV_seqonly_score"] = predictions

# Also include all features in output
output = pd.concat([output, combined], axis=1)

# Save
output.to_csv(output_file, sep="\t", index=False)
sys.stderr.write(f"[scoring_seqonly] Scores written to {output_file}\n")
