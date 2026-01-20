#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Preprocess sequence-only input for SegmentNT.

Input: TSV with columns REF, ALT, [TYPE], [ID]
- REF = reference sequence (required)
- ALT = alternate/variant sequence (required)
- TYPE = SV type (optional, default: "SV")
- ID = variant identifier (optional, default: "var_N")

Processing:
- First 96bp treated as upstream flank
- Last 96bp treated as downstream flank
- Middle part: if > 1008bp, shrink to first 498bp + 12 N's + last 498bp
- N's transformed: groups of 6 N's become 1 N, remainder become A's

Output: Two TSV files for SegmentNT (sample/ALT and reference/REF)
"""

import sys


def transform_genomic_sequence(sequence):
    """
    Transform N runs: groups of 6 N's become 1 N, remainder become A's.
    This matches the tokenization behavior of SegmentNT.
    """
    transformed = []
    i = 0
    while i < len(sequence):
        if sequence[i] == 'N':
            start = i
            while i < len(sequence) and sequence[i] == 'N':
                i += 1
            n_count = i - start
            transformed.append('N' * (n_count // 6))
            transformed.append('A' * (n_count % 6))
        else:
            transformed.append(sequence[i])
            i += 1
    return ''.join(transformed)


def shrink_middle(seq, flank_size=96, max_middle=1008, n_insert=12):
    """
    Process sequence: keep first and last flank_size bp as flanks,
    shrink middle if > max_middle bp.

    For sequences > (2 * flank_size + max_middle):
    - Keep first flank_size as upstream
    - Middle: if > max_middle, take first (max_middle-n_insert)/2 + N*n_insert + last (max_middle-n_insert)/2
    - Keep last flank_size as downstream
    """
    if len(seq) <= 2 * flank_size:
        # Sequence is too short to have distinct flanks and middle
        # Just return as-is
        return seq

    upstream = seq[:flank_size]
    downstream = seq[-flank_size:]
    middle = seq[flank_size:-flank_size]

    if len(middle) > max_middle:
        half = (max_middle - n_insert) // 2  # 498
        middle = middle[:half] + 'N' * n_insert + middle[-half:]

    return upstream + middle + downstream


def adjusted_length(seq: str) -> int:
    """
    Calculate adjusted length where N counts as 6 (matching SegmentNT tokenization).
    This matches the behavior in SB_features_df_generator.py.
    """
    return sum(6 if c == 'N' else 1 for c in seq)


def process_sequence(seq):
    """Apply shrinking and N-transformation to a sequence."""
    seq = seq.upper().strip()
    if not seq:
        # Return minimal valid sequence if empty
        return 'N' * 192  # 96 + 96 for flanks

    # Shrink middle if needed
    processed = shrink_middle(seq)

    # Transform N runs
    if 'N' in processed:
        processed = transform_genomic_sequence(processed)

    return processed


def main():
    if len(sys.argv) != 4:
        sys.stderr.write(
            "Usage: preprocess_sequences.py input.tsv output_alt.bed output_ref.bed\n"
            "Input: TSV with columns REF, ALT, [TYPE], [ID]\n"
        )
        sys.exit(2)

    input_file = sys.argv[1]
    output_alt = sys.argv[2]
    output_ref = sys.argv[3]

    with open(input_file, 'r') as fin, \
         open(output_alt, 'w') as f_alt, \
         open(output_ref, 'w') as f_ref:

        for line_idx, line in enumerate(fin):
            line = line.rstrip('\n')
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 2:
                sys.stderr.write(
                    f"[preprocess_sequences] Skipping line {line_idx}: "
                    f"expected at least 2 columns (REF, ALT), got {len(parts)}\n"
                )
                continue

            # Required columns
            ref_seq = parts[0].upper().strip()
            alt_seq = parts[1].upper().strip()

            # Optional columns with defaults
            sv_type = parts[2] if len(parts) > 2 and parts[2].strip() else "SV"
            id_ = parts[3] if len(parts) > 3 and parts[3].strip() else f"var_{line_idx}"

            # Validate flanking regions match (first and last 96bp)
            flank_size = 96
            ref_upstream = ref_seq[:flank_size] if len(ref_seq) >= flank_size else ref_seq
            alt_upstream = alt_seq[:flank_size] if len(alt_seq) >= flank_size else alt_seq
            ref_downstream = ref_seq[-flank_size:] if len(ref_seq) >= flank_size else ref_seq
            alt_downstream = alt_seq[-flank_size:] if len(alt_seq) >= flank_size else alt_seq

            if ref_upstream != alt_upstream:
                sys.stderr.write(
                    f"[preprocess_sequences] Error line {line_idx} (ID: {id_}): "
                    f"first {flank_size}bp of REF and ALT do not match.\n"
                    f"  REF: {ref_upstream[:50]}{'...' if len(ref_upstream) > 50 else ''}\n"
                    f"  ALT: {alt_upstream[:50]}{'...' if len(alt_upstream) > 50 else ''}\n"
                )
                sys.exit(1)

            if ref_downstream != alt_downstream:
                sys.stderr.write(
                    f"[preprocess_sequences] Error line {line_idx} (ID: {id_}): "
                    f"last {flank_size}bp of REF and ALT do not match.\n"
                    f"  REF: {ref_downstream[:50]}{'...' if len(ref_downstream) > 50 else ''}\n"
                    f"  ALT: {alt_downstream[:50]}{'...' if len(alt_downstream) > 50 else ''}\n"
                )
                sys.exit(1)

            # Process sequences
            processed_alt = process_sequence(alt_seq)
            processed_ref = process_sequence(ref_seq)

            # Output format: ID, 0, len, TYPE, SEQ
            len_alt = len(processed_alt)
            len_ref = len(processed_ref)

            f_alt.write(f"{id_}\t0\t{len_alt}\t{sv_type}\t{processed_alt}\n")
            f_ref.write(f"{id_}\t0\t{len_ref}\t{sv_type}\t{processed_ref}\n")

    sys.stderr.write(f"[preprocess_sequences] Done. Output written to {output_alt} and {output_ref}\n")


if __name__ == "__main__":
    main()
