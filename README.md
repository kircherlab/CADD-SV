# CADD-SV

CADD-SV is a command-line tool for scoring structural variants (SVs). The
`caddsv` command wraps the packaged Snakemake workflow, prepares inputs in the
layout expected by the workflow, runs scoring, and copies the final score tables
to a predictable output directory.

## Quick Start

Install CADD-SV in an isolated environment:

```bash
conda create -n caddsv python=3.12 pip
conda activate caddsv
git clone https://github.com/kircherlab/CADD-SV.git
cd CADD-SV
git checkout dev
pip install .
```

Download the annotation bundle and SegmentNT model files:

```bash
caddsv get annotations \
  --annotations-dir /data/caddsv/annotations \
  --with-segmentnt
```

Score a BED file:

```bash
caddsv run sample.bed \
  --annotations-dir /data/caddsv/annotations \
  --output-dir /data/caddsv/runs/sample \
  --threads 8
```

Final scores are written to:

```text
/data/caddsv/runs/sample/scored/sample_score.tsv
```

Run the sequence-resolved model with SegmentNT-derived features:

```bash
caddsv run sample.bed \
  --seqresolved \
  --annotations-dir /data/caddsv/annotations \
  --output-dir /data/caddsv/runs/sample_seqresolved \
  --threads 8
```

SegmentNT is practical for normal use with a GPU. CPU execution can work for
very small tests, but it is slow.

## Installation

CADD-SV itself is installed with `pip install .` from a source checkout. Use the
quick-start commands above for a fresh install.

The command-line wrapper needs Python and the Python dependencies from
`pyproject.toml`. The pip package installs the CLI and bundled workflow files;
full scoring also requires conda at runtime so Snakemake can create the
workflow environments defined under `caddsv/workflow/envs/`.

The first scoring run may spend time creating those environments under:

```text
${XDG_CACHE_HOME:-$HOME/.cache}/caddsv/snakemake-conda/
```

Later runs reuse the same environments.

To choose a cluster scratch or shared environment location, use either:

```bash
caddsv run sample.bed --conda-prefix /scratch/$USER/caddsv-conda
```

or:

```bash
export CADD_SV_CONDA_PREFIX=/scratch/$USER/caddsv-conda
caddsv run sample.bed
```

Old local `caddsv/.snakemake-envs/` directories from source checkouts are
legacy runtime files and can be removed manually when no run is using them.

## Recommended Layout

Use explicit paths for annotations and outputs when running from different
working directories:

```text
/data/caddsv/
  annotations/
    CADD/
    ucsc/
    segment_nt/
  runs/
    sample/
    sample_seqresolved/
```

If `--annotations-dir` is not provided, CADD-SV uses `./annotations` relative to
the directory where the command is executed. If `--output-dir` is not provided,
CADD-SV uses `./caddsv_results`.

## Downloading Data

### Annotations

Download the v2.0 annotation dependency archive:

```bash
caddsv get annotations
```

By default, files are extracted to `./annotations`. To use a stable shared path:

```bash
caddsv get annotations --annotations-dir /data/caddsv/annotations
```

The archive is downloaded from:

```text
https://kircherlab.bihealth.org/download/CADD-SV/v2.0/dependencies.tar.gz
```

Use the same annotation directory for download and scoring:

```bash
caddsv get annotations --annotations-dir /data/caddsv/annotations
caddsv run sample.bed --annotations-dir /data/caddsv/annotations
```

For shared filesystems, containers, or job schedulers, make the whole annotation
directory available at runtime.

### SegmentNT

SegmentNT model files are used by `--seqresolved` and `--seqonly`. Download them
once into the annotation directory:

```bash
caddsv get segmentnt --annotations-dir /data/caddsv/annotations
```

This creates:

```text
/data/caddsv/annotations/segment_nt/
```

You can download annotations and SegmentNT in one command:

```bash
caddsv get annotations \
  --annotations-dir /data/caddsv/annotations \
  --with-segmentnt
```

When `<annotations-dir>/segment_nt` exists, CADD-SV automatically points
SegmentNT to that local directory and prevents runtime model downloads.

To replace an existing local SegmentNT directory:

```bash
caddsv get segmentnt \
  --annotations-dir /data/caddsv/annotations \
  --force-segmentnt
```

To use a model directory outside the annotation bundle:

```bash
SEGMENTNT_MODEL=/models/segment_nt \
caddsv run sample.bed --seqresolved --annotations-dir /data/caddsv/annotations
```

For fully offline runs, use a local model directory and set:

```bash
HF_HUB_OFFLINE=1
TRANSFORMERS_OFFLINE=1
SEGMENTNT_LOCAL_FILES_ONLY=1
```

SegmentNT is downloaded from `InstaDeepAI/segment_nt` on Hugging Face and is
licensed separately from CADD-SV under CC BY-NC-SA 4.0. Commercial use and
redistribution are subject to the upstream SegmentNT license.

## Running CADD-SV

### Coordinate-Based Scoring

Score one BED file:

```bash
caddsv run sample.bed \
  --annotations-dir /data/caddsv/annotations
```

Score multiple BED files in one invocation:

```bash
caddsv run sample1.bed sample2.bed \
  --annotations-dir /data/caddsv/annotations \
  --output-dir batch_results
```

Final output files are copied to:

```text
<output-dir>/scored/
```

### Sequence-Resolved Scoring

`--seqresolved` adds SegmentNT-derived features to coordinate-based scoring:

```bash
caddsv run sample.bed \
  --seqresolved \
  --annotations-dir /data/caddsv/annotations \
  --output-dir sample_seqresolved
```

This mode still needs coordinate annotations, so keep `--annotations-dir` pointed
at the full annotation bundle. It also needs SegmentNT either in
`<annotations-dir>/segment_nt` or in the directory specified by
`SEGMENTNT_MODEL`.

### Sequence-Only Scoring

`--seqonly` scores variants from paired reference and alternate sequences rather
than genomic coordinates:

```bash
caddsv run sequences.tsv \
  --seqonly \
  --annotations-dir /data/caddsv/annotations \
  --output-dir seqonly_results
```

Sequence-only mode needs SegmentNT, but it does not use the coordinate
annotation tracks. If the model is not under the default `./annotations`
directory, pass the annotation directory that contains `segment_nt` or set
`SEGMENTNT_MODEL`.

Final sequence-only output:

```text
seqonly_results/scored/sequences_seqonly_score.tsv
```

### Reusing Prepared Inputs

When a BED file is passed to `caddsv run`, the wrapper writes a normalized copy
to:

```text
<output-dir>/input/id_<dataset>.bed
```

You can later pass the dataset name instead of the original file:

```bash
caddsv run sample \
  --output-dir caddsv_results \
  --annotations-dir /data/caddsv/annotations
```

This expects:

```text
caddsv_results/input/id_sample.bed
```

For sequence-only mode, the expected file is:

```text
caddsv_results/input/id_sample.tsv
```

## Inputs

### BED Input

Coordinate-based modes use uncompressed `.bed` files with at least four
tab-separated columns:

```text
chrom    start    end    type
```

Supported SV types:

```text
DEL DUP INS INV
```

Example:

```text
chr1    1000000    1001000    DEL
chr1    2000000    2000500    DUP
chr2    3000000    3000000    INS
chr3    4000000    4001200    INV
```

The CLI preprocesses BED input before running Snakemake:

- Adds `chr` to chromosome names when missing.
- Keeps only `chr1` through `chr22`, `chrX`, and `chrY`.
- Keeps only `DEL`, `DUP`, `INS`, and `INV`.
- Skips rows with fewer than four columns.
- Sorts by chromosome and start coordinate.
- Writes normalized input to `<output-dir>/input/id_<dataset>.bed`.

The dataset name is derived from the filename. For `sample.bed`, the dataset is
`sample`. For `id_sample.bed`, the `id_` prefix is stripped and the dataset is
still `sample`.

The wrapper recognizes BED input by the `.bed` suffix. Compressed `.bed.gz`
files are not auto-preprocessed by the wrapper; decompress them first or prepare
`<output-dir>/input/id_<dataset>.bed` manually and run by dataset name.

If the normalized target already exists with different content, CADD-SV asks
before overwriting it. Use a different `--output-dir` when comparing inputs with
the same filename stem.

### Sequence-Only TSV Input

`--seqonly` requires `.tsv` input with positional columns. Each row must have at
least `REF` and `ALT`; `TYPE` and `ID` are optional. Do not include a header row
unless it is an actual sequence record.

```text
REF    ALT    TYPE    ID
```

Column behavior:

| Column | Required | Default |
| --- | --- | --- |
| `REF` | Yes | None |
| `ALT` | Yes | None |
| `TYPE` | No | `SV` |
| `ID` | No | `var_<line_number>` |

Sequence-only preprocessing:

- Converts REF and ALT to uppercase.
- Requires the first 96 bp of REF and ALT to match.
- Requires the last 96 bp of REF and ALT to match.
- Shrinks long middle sequence while preserving flanks.
- Transforms `N` runs to match SegmentNT tokenization.

## Outputs

For `sample.bed` and the default output directory:

```text
caddsv_results/
  input/
    id_sample.bed
  beds/
    sample/
      output/
        samplebed_score100.bed
      ...
  scored/
    sample_score.tsv
```

For `sequences.tsv --seqonly`:

```text
caddsv_results/
  input/
    id_sequences.tsv
  beds/
    sequences/
      output/
        sequences_seqonly_score.tsv
      ...
  scored/
    sequences_seqonly_score.tsv
```

The `scored/` directory is the stable user-facing output location. The `beds/`
directory contains Snakemake intermediates and the workflow's native output.

Main score columns:

| Mode | Main score columns |
| --- | --- |
| Coordinate scoring | `CADD-SV_PHRED`, `CADD-SV_score` |
| Sequence-resolved scoring | `CADD-SV_PHRED`, `CADD-SV_score`, `CADD-SV-SR_PHRED`, `CADD-SV-SR_score` |
| Sequence-only scoring | `CADD-SV_seqonly_PHRED`, `CADD-SV_seqonly_score` |

The output also retains annotation and model feature columns for downstream
inspection.

## Options

### `caddsv get`

```bash
caddsv get annotations [--annotations-dir PATH] [--with-segmentnt] [--force-segmentnt]
caddsv get segmentnt   [--annotations-dir PATH] [--force-segmentnt] [--segmentnt-repo REPO]
```

| Option | Meaning |
| --- | --- |
| `--annotations-dir PATH` | Destination annotation directory. Default: `./annotations`. |
| `--with-segmentnt` | Also download SegmentNT into `<annotations-dir>/segment_nt` while downloading annotations. |
| `--force-segmentnt` | Replace files in an existing local SegmentNT directory before downloading. |
| `--segmentnt-repo REPO` | Hugging Face repository for SegmentNT files. Default: `InstaDeepAI/segment_nt`. |

### `caddsv run`

```bash
caddsv run INPUT [INPUT ...] [OPTIONS]
```

| Option | Meaning |
| --- | --- |
| `--threads`, `-j` | Maximum Snakemake jobs. Default: `4`. |
| `--annotations-dir PATH` | Annotation directory. Default: `./annotations`. |
| `--output-dir`, `-o PATH` | Results directory. Default: `./caddsv_results`. |
| `--conda-prefix PATH` | Snakemake conda environment directory. Default: `CADD_SV_CONDA_PREFIX` or `${XDG_CACHE_HOME:-$HOME/.cache}/caddsv/snakemake-conda`. |
| `--config`, `-c PATH` | Optional YAML configuration file. Default: packaged `caddsv/config.yml`. |
| `--seqresolved` | Add SegmentNT-derived features to coordinate-based scoring. |
| `--seqonly` | Run sequence-only scoring from REF/ALT TSV input. |
| `--force` | Pass `--forceall` to Snakemake. |
| `--unlock` | Unlock a locked Snakemake output directory. |
| `--check-time` | Write a `caddsv_run_<timestamp>.log` resource summary. |

## Runtime Notes

- First runs are slower because Snakemake creates conda environments.
- Snakemake conda environments are stored in the user cache by default; set
  `CADD_SV_CONDA_PREFIX` or `--conda-prefix` to use scratch/shared storage.
- `--threads` controls Snakemake cores, but some steps are I/O-bound.
- SegmentNT is much faster on GPU than CPU.
- Keep annotations and outputs on fast local storage when possible.
- Use the same `--output-dir` to resume or reuse work from an interrupted run.
- Use a new `--output-dir` when comparing inputs with the same filename stem.

To remove cached Snakemake conda environments:

```bash
rm -rf "${XDG_CACHE_HOME:-$HOME/.cache}/caddsv/snakemake-conda"
```

Use `--check-time` to record a small resource summary:

```bash
caddsv run sample.bed \
  --annotations-dir /data/caddsv/annotations \
  --check-time
```

This writes:

```text
caddsv_run_<YYYYMMDD_HHMMSS>.log
```

The log includes the Snakemake command, return code, wall time, CPU time, CPU
utilization, and maximum RSS reported by the operating system.

## Configuration

Most users should prefer CLI flags over editing config files. Use `--config`
only when you need to provide an alternate Snakemake YAML configuration:

```bash
caddsv run sample.bed --config custom.yml
```

The packaged default config is `caddsv/config.yml`. Runtime values such as the
dataset name, annotation directory, output directory, and SegmentNT model path
are normally supplied by the CLI.

## Troubleshooting

### Missing Annotations

Download annotations and run with the same path:

```bash
caddsv get annotations --annotations-dir /data/caddsv/annotations
caddsv run sample.bed --annotations-dir /data/caddsv/annotations
```

### SegmentNT Downloads at Runtime

Download SegmentNT locally:

```bash
caddsv get segmentnt --annotations-dir /data/caddsv/annotations
```

Then run with that annotation directory:

```bash
caddsv run sample.bed --seqresolved --annotations-dir /data/caddsv/annotations
```

Or set the model directory explicitly:

```bash
SEGMENTNT_MODEL=/data/caddsv/annotations/segment_nt \
caddsv run sample.bed --seqresolved --annotations-dir /data/caddsv/annotations
```

### Locked Snakemake Directory

Unlock the output directory:

```bash
caddsv run sample.bed --unlock --output-dir caddsv_results
```

Then rerun the original command.

### Existing Input Prompt

If `<output-dir>/input/id_<dataset>.bed` or `.tsv` exists with different content,
CADD-SV asks before overwriting. Choose a new `--output-dir` to avoid prompts
when comparing multiple inputs with the same dataset name.

### Slow First Run

Common causes:

- Snakemake is creating conda environments.
- SegmentNT or PyTorch dependencies are being installed.
- SegmentNT is running on CPU.
- Annotation files are being read from slow storage.

## Minimal Smoke Test

Create a tiny BED file:

```bash
cat > one_del.bed <<'EOF'
chr1	1000000	1000010	DEL
EOF
```

Run coordinate-based scoring:

```bash
caddsv run one_del.bed \
  --annotations-dir /data/caddsv/annotations \
  --output-dir test_run \
  --threads 1 \
  --check-time
```

Expected final output:

```text
test_run/scored/one_del_score.tsv
```

Run sequence-resolved scoring if SegmentNT is downloaded:

```bash
caddsv run one_del.bed \
  --seqresolved \
  --annotations-dir /data/caddsv/annotations \
  --output-dir test_seqresolved \
  --threads 1 \
  --check-time
```

Expected final output:

```text
test_seqresolved/scored/one_del_score.tsv
```
