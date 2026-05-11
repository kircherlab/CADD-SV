#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
from pathlib import Path

import numpy as np
import pandas as pd

MODEL_DIR = Path(__file__).resolve().parent.parent / "models"
MODEL_PATH = MODEL_DIR / "seqonlymodel.joblib"
META_COLS = ["chr", "start", "end", "type"]


def load_model():
    import joblib

    return joblib.load(MODEL_PATH)


def score_seqonly(sb_file, sbref_file, db_file, output_file, model=None):
    model = model or load_model()

    sb_df = pd.read_table(sb_file, low_memory=False)
    sbref_df = pd.read_table(sbref_file, low_memory=False)
    db_df = pd.read_table(db_file, low_memory=False)

    meta = sb_df[META_COLS]
    sb_features = sb_df.drop(columns=META_COLS)
    sbref_features = sbref_df.drop(columns=META_COLS)
    db_features = db_df.drop(columns=META_COLS)
    sbref_features.columns = [f"{col}_alt" for col in sbref_features.columns]

    features = pd.concat([sb_features, sbref_features, db_features], axis=1)
    features["copies"] = np.where(meta["type"].values == "DEL", -1, 1)

    predictions = model.predict_proba(features[model.feature_names_in_])[:, 1]

    output_parts = []
    ids = meta["chr"].replace("", np.nan).reset_index(drop=True)
    if ids.notna().any():
        output_parts.append(ids.rename("id"))
    output_parts.append(meta["type"].reset_index(drop=True))

    output = pd.concat(output_parts, axis=1)
    output["CADD-SV_seqonly_score"] = predictions
    output = pd.concat([output, features], axis=1)
    output.to_csv(output_file, sep="\t", index=False)
    sys.stderr.write(f"[scoring_seqonly] Scores written to {output_file}\n")


def main():
    if len(sys.argv) != 5:
        sys.stderr.write("Usage: scoring_seqonly.py SB.bed SBref.bed DB.bed output.tsv\n")
        sys.exit(2)

    score_seqonly(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])


if __name__ == "__main__":
    main()
