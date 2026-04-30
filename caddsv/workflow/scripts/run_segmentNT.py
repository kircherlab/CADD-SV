#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
SegmentNT inference with fixed windows of 200 DNA tokens (excluding specials).

Input : TSV file; DNA sequence is in a selectable 1-based column.
Output: HDF5; one group per input line with:
    - probs: float32 array [seq_len_nt, num_features]
    - features: feature names (order matches probs' last dim)
    - attrs: seq_len_nt, num_features, tokens_per_window_excl_specials=200,
             specials_added, window_nt

Progress: prints updates like [current/total] as it processes lines.
"""

import argparse
import os
from pathlib import Path
import sys
import h5py
import numpy as np
import torch
from transformers import AutoTokenizer, AutoModel

DEFAULT_SEGMENTNT_MODEL = "InstaDeepAI/segment_nt"


def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("input_tsv", type=str)
    p.add_argument("output_h5", type=str)
    p.add_argument("--seq-col", type=int, default=5,
                   help="1-based column index of the DNA sequence in the TSV (default: 5)")
    p.add_argument("--device", type=str, default=None, choices=[None, "cpu", "cuda"],
                   help="Force device; default auto-detect.")
    return p.parse_args()


def env_flag(name: str, default: bool = False) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    return value.strip().lower() not in {"0", "false", "no", "off"}


def ensure_pad_token(tokenizer, model):
    """Add a pad token if missing (some tokenizers do not define one)."""
    if getattr(tokenizer, "pad_token", None) is None:
        tokenizer.add_special_tokens({"pad_token": "[PAD]"})
        try:
            model.resize_token_embeddings(len(tokenizer))
        except Exception:
            pass


def count_specials(tokenizer) -> int:
    """Return the number of special tokens added for a single sequence."""
    try:
        return int(tokenizer.num_special_tokens_to_add(pair=False))
    except Exception:
        return 1  # conservative fallback


def tokenize_fixed_len(tokenizer, seq: str, dna_tokens_excl_specials: int, specials: int):
    """
    Produce exactly (dna_tokens_excl_specials + specials) tokens using padding and truncation.
    """
    max_length = dna_tokens_excl_specials + specials
    return tokenizer(
        seq,
        return_tensors="pt",
        padding="max_length",
        truncation=True,
        max_length=max_length,
    )


def feature_names(model, logits=None):
    names = getattr(model.config, "features", None)
    if names:
        try:
            return list(names)
        except Exception:
            pass
    if logits is not None:
        num_features = int(logits.shape[-2])  # [B, L_nt, F, 2]
        return [f"feature_{i}" for i in range(num_features)]
    return None


def main():
    args = parse_args()
    seq_col0 = args.seq_col - 1

    model_source = os.environ.get("SEGMENTNT_MODEL", DEFAULT_SEGMENTNT_MODEL)
    local_files_only = (
        env_flag("SEGMENTNT_LOCAL_FILES_ONLY")
        or env_flag("HF_HUB_OFFLINE")
        or env_flag("TRANSFORMERS_OFFLINE")
        or Path(model_source).exists()
    )

    tokenizer = AutoTokenizer.from_pretrained(
        model_source,
        trust_remote_code=True,
        local_files_only=local_files_only,
    )
    model = AutoModel.from_pretrained(
        model_source,
        trust_remote_code=True,
        local_files_only=local_files_only,
    )

    # Device
    device = torch.device(args.device or ("cuda" if torch.cuda.is_available() else "cpu"))
    model.to(device).eval()
    ensure_pad_token(tokenizer, model)

    # Strict window: 200 DNA tokens excluding specials. Must be divisible by 4.
    DNA_TOKENS_EXCL_SPECIALS = 200
    assert DNA_TOKENS_EXCL_SPECIALS % 4 == 0
    specials = count_specials(tokenizer)
    window_tokens_total = DNA_TOKENS_EXCL_SPECIALS + specials

    # SegmentNT uses non-overlapping 6-mers -> ~6 nt per DNA token
    window_nt = DNA_TOKENS_EXCL_SPECIALS * 6

    # Pre-count total lines we will actually process (non-empty and with the sequence column present)
    with open(args.input_tsv, "r") as fin:
        total_to_process = 0
        for _line in fin:
            if not _line.strip():
                continue
            cols = _line.rstrip("\n").split("\t")
            if seq_col0 < len(cols):
                total_to_process += 1

    processed = 0

    with open(args.input_tsv, "r") as fin, h5py.File(args.output_h5, "w") as h5:
        seen_feature_names = None

        for file_line_idx, line in enumerate(fin):
            line = line.rstrip("\n")
            if not line:
                continue

            cols = line.split("\t")
            if seq_col0 >= len(cols):
                continue

            # Progress update
            processed += 1
            sys.stderr.write(f"\r[{processed}/{total_to_process}]")
            sys.stderr.flush()

            seq = cols[seq_col0].strip().upper()
            if not seq:
                # ensure at least one window
                seq = "N" * window_nt

            seq_len_nt = len(seq)
            probs_chunks = []

            # Tile sequence into 200-token windows (~1200 nt)
            for start in range(0, seq_len_nt, window_nt):
                chunk = seq[start:start + window_nt]
                real_len_nt = len(chunk)

                toks = tokenize_fixed_len(tokenizer, chunk, DNA_TOKENS_EXCL_SPECIALS, specials)
                input_ids = toks["input_ids"].to(device)
                attention_mask = toks["attention_mask"].to(device)

                with torch.no_grad():
                    out = model(input_ids=input_ids, attention_mask=attention_mask)

                # logits: [B, L_nt_window, F, 2]
                logits = out.logits

                if seen_feature_names is None:
                    seen_feature_names = feature_names(model, logits)

                # Convert to probabilities; keep class index 1
                probs = torch.softmax(logits, dim=-1)[..., 1]   # [1, L_nt_window, F]
                probs = probs[:, :real_len_nt, :]               # drop positions from padding
                probs_chunks.append(probs.squeeze(0).cpu().numpy().astype(np.float32))

                if device.type == "cuda":
                    torch.cuda.empty_cache()

            if probs_chunks:
                probs_full = np.concatenate(probs_chunks, axis=0)  # [seq_len_nt, F]
            else:
                f = len(seen_feature_names or [])
                probs_full = np.zeros((0, f), dtype=np.float32)

            # Save one group per line (use the source file line index for traceability)
            grp = h5.create_group(f"group_{file_line_idx}")
            grp.create_dataset("probs", data=probs_full, compression="gzip")

            # Save feature names (UTF-8 strings)
            names = seen_feature_names or [f"feature_{i}" for i in range(probs_full.shape[1])]
            str_dt = h5py.string_dtype(encoding="utf-8")
            grp.create_dataset("features", data=np.array(names, dtype=object), dtype=str_dt)

            # Metadata
            grp.attrs["seq_len_nt"] = int(probs_full.shape[0])
            grp.attrs["num_features"] = int(probs_full.shape[1]) if probs_full.size else 0
            grp.attrs["tokens_per_window_excl_specials"] = int(DNA_TOKENS_EXCL_SPECIALS)
            grp.attrs["specials_added"] = int(specials)
            grp.attrs["window_nt"] = int(window_nt)

    # Finish the progress line with a newline
    sys.stderr.write("\n")
    sys.stderr.flush()


if __name__ == "__main__":
    main()
