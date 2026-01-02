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

This command will:

1. Download a pre-packaged archive of annotation dependencies.
2. Unpack the archive into the working directory.

This step only needs to be performed once per installation.

---

## Running CADD-SV

### Basic Usage

To score a BED file containing structural variants:

```bash
caddsv run file.bed
```

---

### Sequence Model (SegmentNT)

To run the model with integration of SegmentNT-derived annotations, enable the sequence-based model:

```bash
caddsv run file.bed --sequence_model
```

This flag propagates directly into the Snakemake configuration as `sequence_model=True`.

---

### Threads

Control the maximum number of parallel jobs:

```bash
caddsv run file.bed --threads 8
```

(Default: 4)

---

### Configuration File

An optional YAML configuration file can be supplied:

```bash
caddsv run file.bed --config custom_config.yml
```

If not provided, the default packaged `config.yml` is used.

---

### Forcing or Unlocking Runs

- Force rerun of all Snakemake rules:
  ```bash
  caddsv run file.bed --force
  ```

- Unlock a previously locked Snakemake working directory:
  ```bash
  caddsv run file.bed --unlock
  ```

---

## Input Handling and Preprocessing

### BED File Requirements

When a BED file is supplied, the CLI performs preprocessing before invoking Snakemake:

- Ensures chromosome names are prefixed with `chr`
- Filters to allowed chromosomes:
  - `chr1`–`chr22`, `chrX`, `chrY`
- Filters to allowed SV types (4th column):
  - `DEL`, `DUP`, `INS`, `INV`
- Requires at least four columns per row
- Removes invalid or malformed entries
- Sorts the BED file by chromosome and start coordinate

Processed BED files are written to:

```text
input/id_<dataset>.bed
```

---

## Output

Final scores are written to the `scored/` directory:

```text
scored/<dataset>_score.tsv
```
