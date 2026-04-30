# CADD-SV

CADD-SV is a Snakemake-based command-line tool for scoring structural variants (SVs).
The provided CLI (`caddsv`) is a lightweight wrapper around an underlying Snakemake workflow, handling input normalization, validation, and execution orchestration.

The primary user-facing functionality is exposed via two commands:

- `caddsv get annotations`
- `caddsv run`

---

## Installation

Clone the repository and install the package using `pip`:

```bash
git clone https://github.com/kircherlab/CADD-SV.git
cd CADD-SV
git checkout dev
pip install .
```

---

## Downloading Annotations

Before running scoring, required annotation dependencies must be downloaded and unpacked.

```bash
caddsv get annotations
```

By default, annotations are extracted to `./annotations`. To specify a custom location:

```bash
caddsv get annotations --annotations-dir /path/to/annotations
```

This step only needs to be performed once per installation.

### Downloading SegmentNT Locally

SegmentNT model files can be downloaded into the annotation directory so
sequence-resolved and sequence-only runs do not fetch the model at runtime:

```bash
caddsv get segmentnt --annotations-dir /path/to/annotations
```

This creates:

```text
/path/to/annotations/segment_nt/
```

You can also download annotations and SegmentNT in one command:

```bash
caddsv get annotations --annotations-dir /path/to/annotations --with-segmentnt
```

When `annotations/segment_nt` exists, CADD-SV automatically runs SegmentNT from
that local folder. To override the model location manually:

```bash
SEGMENTNT_MODEL=/path/to/segment_nt caddsv run file.bed --seqresolved
```

The SegmentNT files are downloaded from `InstaDeepAI/segment_nt` on Hugging Face
and are licensed separately from CADD-SV under CC BY-NC-SA 4.0.

---

## Running CADD-SV

### Basic Usage

To score a BED file containing structural variants:

```bash
caddsv run file.bed
```

Multiple BED files can be scored in a single invocation:

```bash
caddsv run file1.bed file2.bed
```

Results are written to `caddsv_results/scored/` by default.

---

### Sequence-Resolved Model (SegmentNT)

To run the model with integration of SegmentNT-derived annotations (requires GPU):

```bash
caddsv run file.bed --seqresolved
```

---

### Sequence-Only Mode

To score variants using only SegmentNT-derived features, without coordinate-based annotations:

```bash
caddsv run sequences.tsv --seqonly
```

Input must be a TSV file with columns: `REF`, `ALT`, `[TYPE]`, `[ID]`.

Scores are written to `caddsv_results/scored/<dataset>_seqonly_score.tsv`.

---

### Options

| Flag | Description |
|------|-------------|
| `--threads`, `-j` | Max parallel jobs (default: 4) |
| `--config`, `-c` | Optional YAML configuration file (default: packaged `config.yml`) |
| `--output-dir`, `-o` | Results directory (default: `./caddsv_results`) |
| `--annotations-dir` | Path to annotation directory (default: `./annotations`) |
| `--force` | Force rerun of all Snakemake rules |
| `--unlock` | Unlock a previously locked Snakemake working directory |
| `--check-time` | Track time and resource usage; log to a `.log` file |

---

## Input Handling and Preprocessing

### BED File Requirements

When a BED file is supplied, the CLI performs preprocessing before invoking Snakemake:

- Ensures chromosome names are prefixed with `chr`
- Filters to allowed chromosomes: `chr1`--`chr22`, `chrX`, `chrY`
- Filters to allowed SV types (4th column): `DEL`, `DUP`, `INS`, `INV`
- Requires at least four columns per row
- Removes invalid or malformed entries
- Sorts the BED file by chromosome and start coordinate

Processed files are written to:

```text
caddsv_results/input/id_<dataset>.bed
```

If the target file already exists with different content, CADD-SV will prompt before overwriting.

---

## Output

Final scores are written to:

```text
caddsv_results/scored/<dataset>_score.tsv
```

Output columns: `chr | start | end | type | CADD-SV_PHRED | CADD-SV_score | feature1 | feature2 | ...`

When using `--seqresolved`, additional columns `CADD-SV-SR_PHRED` and `CADD-SV-SR_score` are included.

When using `--seqonly`, output is written to `<dataset>_seqonly_score.tsv`.
