# CADD-SV
[![Release](https://img.shields.io/badge/release-v2.0-green)](https://github.com/kircherlab/CADD-SV/releases/tag/v2.0)
[![PyPI version](https://img.shields.io/pypi/v/caddsv.svg)](https://pypi.org/project/caddsv/)
[![Bioconda version](https://img.shields.io/conda/vn/bioconda/caddsv.svg?style=flat)](https://bioconda.github.io/recipes/caddsv/README.html)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://github.com/kircherlab/CADD-SV/blob/main/LICENSE)

CADD-SV is a command-line tool for scoring structural variants (SVs). The
`caddsv` command wraps the packaged Snakemake workflow, prepares input files,
runs scoring, and copies final score tables into a stable output directory.

## Quick Start

Install from a source checkout:

```bash
conda create -n caddsv python=3.12 pip
conda activate caddsv
git clone https://github.com/kircherlab/CADD-SV.git
cd CADD-SV
pip install .
```

Alternatively, install CADD-SV from Bioconda:

```bash
conda install -c bioconda caddsv
```
Installation through PyPI is also available, conda will still be needed to run it

```bash
pip install caddsv
```
Download the annotation bundle:

```bash
caddsv get annotations --annotations-dir /data/caddsv/annotations
```

Score a BED file:

```bash
caddsv run examples/variants.bed \
  --annotations-dir /data/caddsv/annotations \
  --output-dir /data/caddsv/runs/variants \
  --threads 8
```

Final scores are copied to:

```text
/data/caddsv/runs/variants/scored/variants_score.tsv
```

To run SegmentNT-backed modes, download the model files once:

```bash
caddsv get segmentnt --annotations-dir /data/caddsv/annotations
```

You can also download annotations and SegmentNT together:

```bash
caddsv get annotations \
  --annotations-dir /data/caddsv/annotations \
  --with-segmentnt
```

## Installation Notes

CADD-SV installs with `pip install .` from this repository or from Bioconda. The
package includes the CLI and workflow files. By default, full scoring also uses
conda at runtime so Snakemake can create isolated environments for the workflow
rules.

By default, those environments are cached under:

```text
${XDG_CACHE_HOME:-$HOME/.cache}/caddsv/snakemake-conda/
```

Use `--conda-prefix` or `CADD_SV_CONDA_PREFIX` to place them on scratch or
shared storage:

```bash
caddsv run sample.bed --conda-prefix /scratch/$USER/caddsv-conda
```

For a container or another pre-provisioned environment, disable Snakemake's
per-rule conda environments:

```bash
caddsv run sample.bed --no-use-conda
```

In this mode, CADD-SV does not create a conda cache and does not pass
`--use-conda` or `--conda-prefix` to Snakemake. Every executable and Python
package required by the selected workflow rules must already be available in
the parent environment. The files under `caddsv/workflow/envs/` describe those
rule dependencies.

For Singularity/Apptainer deployments, use the versioned OCI image assigned to
each workflow environment:

```bash
caddsv run sample.bed \
  --use-singularity \
  --singularity-prefix /scratch/$USER/caddsv-singularity
```

Images are cached in `${XDG_CACHE_HOME:-$HOME/.cache}/caddsv/snakemake-singularity`
by default. Override the image URIs, including with local SIF paths on an
air-gapped cluster, through the `containers` mapping in a config file. The
package includes the Dockerfile and environment definitions used to build the
images; image binaries are published separately.

For GPU-enabled SegmentNT execution, pass the runtime flag explicitly:

```bash
caddsv run sample.bed --use-singularity --singularity-args=--nv
```

## Data

### Annotations

The annotation bundle is downloaded from:

```text
https://kircherlab.bihealth.org/download/CADD-SV/v2.0/dependencies.tar.gz
```

The default destination is `./annotations`. For reproducible runs, use an
explicit path and pass the same path to `caddsv run`:

```bash
caddsv get annotations --annotations-dir /data/caddsv/annotations
caddsv run sample.bed --annotations-dir /data/caddsv/annotations
```

### SegmentNT

`--seqresolved` and `--seqonly` require SegmentNT model files. The default local
location is:

```text
<annotations-dir>/segment_nt/
```

If the model lives somewhere else, set `SEGMENTNT_MODEL`:

```bash
SEGMENTNT_MODEL=/models/segment_nt \
caddsv run sample.bed --seqresolved --annotations-dir /data/caddsv/annotations
```

For offline runs, point to a local model directory and set:

```bash
HF_HUB_OFFLINE=1
TRANSFORMERS_OFFLINE=1
SEGMENTNT_LOCAL_FILES_ONLY=1
```

SegmentNT is downloaded from `InstaDeepAI/segment_nt` on Hugging Face and is
licensed separately under CC BY-NC-SA 4.0.

## Recommended Layout

Use explicit annotation and output paths when running from different working
directories:

```text
/data/caddsv/
  annotations/
    CADD/
    ucsc/
    segment_nt/
  runs/
    sample/
```

If paths are omitted, CADD-SV uses `./annotations` and `./caddsv_results`
relative to the current working directory.

## Running CADD-SV

### Coordinate-Based Scoring

```bash
caddsv run examples/variants.bed \
  --annotations-dir /data/caddsv/annotations \
  --output-dir sample_results
```

Multiple BED files can be scored in one invocation:

```bash
caddsv run sample1.bed sample2.bed \
  --annotations-dir /data/caddsv/annotations \
  --output-dir batch_results
```

### Sequence-Resolved Scoring

`--seqresolved` adds SegmentNT-derived features to coordinate-based scoring:

```bash
caddsv run sample.bed \
  --seqresolved \
  --annotations-dir /data/caddsv/annotations \
  --output-dir sample_seqresolved
```

This mode needs both the coordinate annotation bundle and SegmentNT model files.
GPU execution is recommended for normal use; CPU execution is mainly practical
for very small tests.

### Sequence-Only Scoring

`--seqonly` scores REF/ALT sequence pairs instead of genomic coordinates:

```bash
caddsv run examples/sequences.tsv \
  --seqonly \
  --annotations-dir /data/caddsv/annotations \
  --output-dir seqonly_results
```

Sequence-only mode needs SegmentNT. It does not use coordinate annotation
tracks, but `--annotations-dir` is still useful when SegmentNT is stored under
`<annotations-dir>/segment_nt`.

### Reusing Prepared Inputs

When a BED file is passed to `caddsv run`, CADD-SV writes a normalized copy to:

```text
<output-dir>/input/id_<dataset>.bed
```

You can later rerun by dataset name:

```bash
caddsv run sample \
  --output-dir caddsv_results \
  --annotations-dir /data/caddsv/annotations
```

For `--seqonly`, the prepared input is `input/id_<dataset>.tsv`.

## Inputs

### BED

Coordinate-based modes use uncompressed `.bed` files with at least four
tab-separated columns:

```text
chrom    start    end    type    [sequence]
```

BED uses a 0-based start and 1-based end coordinate; interval length is
`end - start`. Supported SV types are `DEL`, `DUP`, `INS`, and `INV`. SVs should
be at least 50 bp; for `INS`, this means providing an inserted sequence of at
least 50 bp in the optional fifth column when running `--seqresolved`.

The repository includes a minimal BED example at `examples/variants.bed`:

```text
chr1    999999     1000049    DEL
chr2    2999999    3000000    INS    ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT
```

Before running Snakemake, the CLI adds missing `chr` prefixes, keeps standard
chromosomes (`chr1` through `chr22`, `chrX`, `chrY`), keeps supported SV types,
skips short rows, sorts by chromosome and start, and writes the normalized file
under `<output-dir>/input/`.

Compressed `.bed.gz` files are not auto-preprocessed; decompress them first or
prepare the normalized input manually.

### Sequence-Only TSV

`--seqonly` requires `.tsv` input with positional columns. Each row must include
`REF` and `ALT`; `TYPE` and `ID` are optional. Do not include a header row unless
it is an actual sequence record.

```text
REF    ALT    TYPE    ID
```

| Column | Required | Default |
| --- | --- | --- |
| `REF` | Yes | None |
| `ALT` | Yes | None |
| `TYPE` | No | `SV` |
| `ID` | No | Blank; omitted from final output when absent |

Sequence-only preprocessing uppercases sequences, requires matching 96 bp
flanks, shrinks long middle sequence, and normalizes `N` runs for SegmentNT
tokenization.

The repository includes a headerless sequence-only example with DEL and INS
records at `examples/sequences.tsv`.

## Outputs

For `sample.bed` and the default output directory:

```text
caddsv_results/
  input/id_sample.bed
  beds/sample/
  scored/sample_score.tsv
```

For `sequences.tsv --seqonly`:

```text
caddsv_results/
  input/id_sequences.tsv
  beds/sequences/
  scored/sequences_seqonly_score.tsv
```

The `scored/` directory is the stable user-facing output location. The `beds/`
directory contains Snakemake intermediates and native workflow outputs.

Main score columns:

| Mode | Main score columns |
| --- | --- |
| Coordinate scoring | `CADD-SV_PHRED`, `CADD-SV_score` |
| Sequence-resolved scoring | `CADD-SV_PHRED`, `CADD-SV_score`, `CADD-SV-SR_PHRED`, `CADD-SV-SR_score` |
| Sequence-only scoring | `CADD-SV_seqonly_PHRED`, `CADD-SV_seqonly_score` |

The output also keeps annotation and model feature columns for downstream
inspection.

## Options

### `caddsv get`

```bash
caddsv get annotations [--annotations-dir PATH] [--with-segmentnt] [--force-segmentnt]
caddsv get segmentnt   [--annotations-dir PATH] [--force-segmentnt] [--segmentnt-repo REPO]
```

| Option | Meaning |
| --- | --- |
| `--annotations-dir PATH` | Annotation directory. Default: `./annotations`. |
| `--with-segmentnt` | Also download SegmentNT into `<annotations-dir>/segment_nt`. |
| `--force-segmentnt` | Replace an existing local SegmentNT directory. |
| `--segmentnt-repo REPO` | Hugging Face SegmentNT repository. Default: `InstaDeepAI/segment_nt`. |

### `caddsv run`

```bash
caddsv run INPUT [INPUT ...] [OPTIONS]
```

| Option | Meaning |
| --- | --- |
| `--threads`, `-j` | Maximum Snakemake jobs. Default: `4`. |
| `--annotations-dir PATH` | Annotation directory. Default: `./annotations`. |
| `--output-dir`, `-o PATH` | Results directory. Default: `./caddsv_results`. |
| `--use-conda` / `--no-use-conda` | Enable or disable Snakemake conda environments. Enabled by default. |
| `--conda-prefix PATH` | Snakemake conda environment directory. |
| `--use-singularity` / `--use-apptainer` | Run rules in their versioned OCI/SIF containers. |
| `--singularity-prefix PATH` | Singularity/Apptainer image cache directory. |
| `--singularity-args TEXT` | Extra runtime arguments, such as `--nv` for GPUs. |
| `--config`, `-c PATH` | Alternate Snakemake YAML configuration. |
| `--seqresolved` | Add SegmentNT-derived features to coordinate-based scoring. |
| `--seqonly` | Run sequence-only scoring from REF/ALT TSV input. |
| `--force` | Pass `--forceall` to Snakemake. |
| `--unlock` | Unlock a locked Snakemake output directory. |
| `--check-time` | Write a small resource summary log. |

## Runtime Notes

- First runs are slower because Snakemake creates or pulls software environments.
- In containerized runs, use `--no-use-conda` for a prebuilt parent environment or `--use-singularity` for per-rule images.
- Use the same `--output-dir` to resume or reuse work from an interrupted run.
- Use a new `--output-dir` when comparing inputs with the same filename stem.
- `--threads` controls Snakemake cores, but some steps are I/O-bound.
- SegmentNT is much faster on GPU than CPU.
- Keep annotations and outputs on fast local storage when possible.

To remove cached Snakemake conda environments:

```bash
rm -rf "${XDG_CACHE_HOME:-$HOME/.cache}/caddsv/snakemake-conda"
```

To record a resource summary:

```bash
caddsv run sample.bed \
  --annotations-dir /data/caddsv/annotations \
  --check-time
```

This writes `caddsv_run_<YYYYMMDD_HHMMSS>.log` with the Snakemake command,
return code, wall time, CPU time, CPU utilization, and maximum RSS.

## Configuration

Most users should prefer CLI flags over editing config files. Use `--config`
only when you need an alternate Snakemake YAML configuration:

```bash
caddsv run sample.bed --config custom.yml
```

The packaged default config is `caddsv/config.yml`.

## Troubleshooting

### Missing Annotations

Download annotations and run with the same path:

```bash
caddsv get annotations --annotations-dir /data/caddsv/annotations
caddsv run sample.bed --annotations-dir /data/caddsv/annotations
```

### SegmentNT Downloads at Runtime

Download SegmentNT locally, then rerun with the same annotation directory:

```bash
caddsv get segmentnt --annotations-dir /data/caddsv/annotations
caddsv run sample.bed --seqresolved --annotations-dir /data/caddsv/annotations
```

If the model is outside the annotation directory, set `SEGMENTNT_MODEL`.

### Locked Snakemake Directory

```bash
caddsv run sample.bed --unlock --output-dir caddsv_results
```

Then rerun the original command.

### Existing Input Prompt

If `<output-dir>/input/id_<dataset>.bed` or `.tsv` exists with different
content, CADD-SV asks before overwriting. Use a new `--output-dir` to avoid
prompts when comparing inputs with the same dataset name.

### Slow First Run

Common causes are conda environment creation, SegmentNT or PyTorch dependency
setup, CPU-based SegmentNT execution, or annotation files on slow storage.

## Minimal Smoke Test

With annotations already downloaded, run the included BED example:

```bash
caddsv run examples/variants.bed \
  --annotations-dir /data/caddsv/annotations \
  --output-dir test_run \
  --threads 1
```

Expected output:

```text
test_run/scored/variants_score.tsv
```

For `--seqresolved`, download SegmentNT first and run:

```bash
caddsv run examples/variants.bed \
  --seqresolved \
  --annotations-dir /data/caddsv/annotations \
  --output-dir test_seqresolved \
  --threads 1
```

Expected output:

```text
test_seqresolved/scored/variants_score.tsv
```

For `--seqonly`, download SegmentNT first and run:

```bash
caddsv run examples/sequences.tsv \
  --seqonly \
  --annotations-dir /data/caddsv/annotations \
  --output-dir test_seqonly \
  --threads 1
```

Expected output:

```text
test_seqonly/scored/sequences_seqonly_score.tsv
```
