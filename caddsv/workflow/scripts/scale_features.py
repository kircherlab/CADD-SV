"""
Compute per-feature z-scores for scored CADD-SV output.

For each variant, the z-score is computed as (value - mean) / std using
per-SV-type statistics tables (e.g. DEL_stats.tsv).  Features absent from
the stats table are silently dropped from the output.
"""

import pandas as pd
import numpy as np
from pathlib import Path

# Columns that are never features
_META_COLS = {"chr", "start", "end", "type"}


def _is_score_col(col: str) -> bool:
    """Return True for score/PHRED columns that should not be z-scaled."""
    low = col.lower()
    return "score" in low or "phred" in low


def _load_stats(stats_dir: Path, sv_type: str) -> pd.DataFrame:
    """Load {SV_TYPE}_stats.tsv and return it indexed by feature name."""
    path = stats_dir / f"{sv_type}_stats.tsv"
    if not path.exists():
        raise FileNotFoundError(
            f"Stats file not found: {path}\n"
            f"Expected a TSV with columns: feature, mean, std"
        )
    return pd.read_table(path).set_index("feature")


def scale_features(scored_path: str, output_path: str, stats_dir: Path) -> None:
    """Read a scored TSV, compute z-scores per SV type, and write output."""
    df = pd.read_table(scored_path, low_memory=False)

    # Identify feature columns (everything that isn't meta or score/PHRED)
    feature_cols = [
        c for c in df.columns if c not in _META_COLS and not _is_score_col(c)
    ]

    # Cache stats tables per SV type
    stats_cache: dict[str, pd.DataFrame] = {}

    # Prepare output: meta columns + feature columns (z-scored)
    meta = df[["chr", "start", "end", "type"]].copy()
    z_df = pd.DataFrame(index=df.index)

    for sv_type in df["type"].unique():
        if sv_type not in stats_cache:
            stats_cache[sv_type] = _load_stats(stats_dir, sv_type)

        stats = stats_cache[sv_type]
        mask = df["type"] == sv_type

        for col in feature_cols:
            if col not in stats.index:
                continue
            mean = stats.loc[col, "mean"]
            std = stats.loc[col, "std"]
            raw = df.loc[mask, col].astype(float)
            if std == 0:
                z_df.loc[mask, col] = 0.0
            else:
                z_df.loc[mask, col] = (raw - mean) / std

    # Keep only feature columns that were actually z-scored
    scaled_cols = [c for c in feature_cols if c in z_df.columns]
    out = pd.concat([meta, z_df[scaled_cols]], axis=1)
    out.to_csv(output_path, sep="\t", index=False)
