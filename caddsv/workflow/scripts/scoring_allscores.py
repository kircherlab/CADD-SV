#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Scoring script for all-scores mode.

Runs all three models:
- standardmodel.joblib (CADD-SV)
- seqresmodel.joblib (CADD-SV-SR)
- seqonlymodel.joblib (CADD-SV-seqonly)

Usage: scoring_allscores.py input.bed XB.bed SB.bed SBref.bed DB.bed output.bed
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

standard_model_path = MODEL_DIR / "standardmodel.joblib"
seqres_model_path = MODEL_DIR / "seqresmodel.joblib"
seqonly_model_path = MODEL_DIR / "seqonlymodel.joblib"

if len(sys.argv) != 7:
    sys.stderr.write(
        "Usage: scoring_allscores.py input.bed XB.bed SB.bed SBref.bed DB.bed output.bed\n"
    )
    sys.exit(2)

input_file = sys.argv[1]
xb_file = sys.argv[2]
sb_file = sys.argv[3]
sbref_file = sys.argv[4]
db_file = sys.argv[5]
output_file = sys.argv[6]

# Load models
rf = joblib.load(standard_model_path)
rfSR = joblib.load(seqres_model_path)
seqonly_model = joblib.load(seqonly_model_path)

# Read main feature file (XB = SBfinal, contains CB + sequence features)
svdata = pd.read_table(xb_file, low_memory=False)
copies_svdata = pd.read_table(input_file, header=None, low_memory=False)
svdata["copies"] = np.where(copies_svdata.iloc[:, 3].values == "DEL", -1, 1)

# Score with sequence-resolved model (CADD-SV-SR)
X_full = svdata[rfSR.feature_names_in_]
pred_sr = rfSR.predict_proba(X_full)[:, 1]

# Score with standard model (CADD-SV)
pred_standard = rf.predict_proba(X_full[rf.feature_names_in_])[:, 1]

# Read SB/SBref/DB feature files for seqonly scoring
sb_df = pd.read_table(sb_file, low_memory=False)
sbref_df = pd.read_table(sbref_file, low_memory=False)
db_df = pd.read_table(db_file, low_memory=False)

# Get feature columns (everything except the ID columns)
id_cols = ["chr", "start", "end", "type"]
sb_features = sb_df.drop(columns=id_cols, errors='ignore')
sbref_features = sbref_df.drop(columns=id_cols, errors='ignore')
db_features = db_df.drop(columns=id_cols, errors='ignore')

# Rename columns to match seqonly model expectations
sb_features.columns = [f"SB_{c}" for c in sb_features.columns]
sbref_features.columns = [f"SBref_{c}" for c in sbref_features.columns]
db_features.columns = [f"DB_{c}" for c in db_features.columns]

# Combine for seqonly model
combined_seqonly = pd.concat([sb_features, sbref_features, db_features], axis=1)

# Score with seqonly model
try:
    X_seqonly = combined_seqonly[seqonly_model.feature_names_in_]
    pred_seqonly = seqonly_model.predict_proba(X_seqonly)[:, 1]
except KeyError as e:
    missing = set(seqonly_model.feature_names_in_) - set(combined_seqonly.columns)
    sys.stderr.write(f"[scoring_allscores] Warning: Missing features for seqonly model: {missing}\n")
    pred_seqonly = np.full(len(copies_svdata), np.nan)

# Build output
copies_svdata.columns = ["chr", "start", "end", "type"]
copies_svdata["CADD-SV_score"] = pred_standard
copies_svdata["CADD-SV-SR_score"] = pred_sr
copies_svdata["CADD-SV-seqonly_score"] = pred_seqonly

output = pd.concat([copies_svdata, svdata.iloc[:, 3:]], axis=1)
output.to_csv(output_file, sep="\t", index=False)

sys.stderr.write(f"[scoring_allscores] All scores written to {output_file}\n")
