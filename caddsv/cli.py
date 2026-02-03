from pathlib import Path
import subprocess
import sys
from typing import List, Optional
import tempfile
import typer
import os
import shutil
import time
from datetime import datetime
import resource

app = typer.Typer(add_completion=False, help="CADD-SV Snakemake-based scoring tool")

PKG_DIR = Path(__file__).resolve().parent
WORKFLOW_DIR = PKG_DIR / "workflow"
DEFAULT_CONFIG = PKG_DIR / "config.yml"


def format_seconds(value: float) -> str:
    return f"{value:.2f}s"


@app.command()
def get(
        flag: str = typer.Argument(None)
):

    if flag == "annotations":
        typer.echo("Downloading dependencies...")
        subprocess.run(["wget",
            "https://kircherlab.bihealth.org/download/CADD-SV/v2.0/dependencies.tar.gz"],
            check=True)
        typer.echo("Uncompressing dependencies...")
        subprocess.run(["tar",
            "-xf",
            "dependencies.tar.gz"],
            check=True)
        typer.echo("DONE")

@app.command()
def run(
    items: List[str] = typer.Argument(
        ...,
        help=(
            "One or more inputs, each either:\n"
            "  - a BED file (e.g. variants.bed), or\n"
            "  - a TSV file with --sequence-only (REF, ALT, [TYPE], [ID]), or\n"
            "  - a dataset name (e.g. longrange_cherie)\n\n"
            "BED files are automatically mapped to input/id_<name>.bed.\n"
            "For names, CADD-SV expects input/id_<name>.bed to exist."
        ),
    ),
    threads: int = typer.Option(4, "--threads", "-j", help="Max parallel jobs"),
    config: Optional[Path] = typer.Option(
        None, "--config", "-c", help="Optional config file (YAML)"
    ),
    mode: str = typer.Option("scoring", help="Mode; default is scoring"),
    force: bool = typer.Option(
        False, "--force", help="Force rerun of all rules (snakemake --forceall)"
    ),
    sequence_model: bool = typer.Option(False, "--sequence_model", help="Optional to run the model with the integration of SegmentNT derived annotations"),
    sequence_only: bool = typer.Option(
        False, "--sequence-only",
        help="Sequence-only mode: input TSV with columns REF, ALT, [TYPE], [ID]. "
             "Runs only SegmentNT, generates SB/SBref/DB features without coordinate-based annotations."
    ),
    unlock: bool = typer.Option(
        False, "--unlock", help="snakemake --unlock"
    ),
    check_time: bool = typer.Option(
        False,
        "--check-time",
        help="Track time and resource usage; log to a .log file.",
    ),
):
    cfg = config if config is not None else DEFAULT_CONFIG

    # Override mode if sequence_only flag is set
    if sequence_only:
        mode = "seqonly"

    workdir = Path("beds")
    workdir.mkdir(exist_ok=True)

    outdir = Path("scored")
    outdir.mkdir(exist_ok=True)

    input_dir = Path("input")
    input_dir.mkdir(exist_ok=True)

    datasets: List[str] = []

    allowed_chroms = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}
    allowed_svtypes = {"DEL", "DUP", "INS", "INV"}

    for raw in items:
        p = Path(raw)

        # Handle sequence-only mode (TSV input)
        if sequence_only:
            if p.suffix.lower() != ".tsv":
                raise typer.BadParameter(
                    f"Sequence-only mode requires .tsv files, got '{p}'"
                )
            if not p.exists():
                raise typer.BadParameter(f"TSV file '{p}' does not exist.")

            tsv_path = p.resolve()
            stem = tsv_path.stem
            if stem.startswith("id_"):
                name = stem[3:]
            else:
                name = stem

            target = input_dir / f"id_{name}.tsv"

            # For TSV, just copy if different (no sorting needed)
            if target.exists():
                with open(tsv_path) as new_f:
                    new_content = new_f.read()
                with open(target) as existing_f:
                    existing_content = existing_f.read()

                if new_content == existing_content:
                    typer.echo(f"Using existing file: {target}")
                else:
                    typer.echo(f"File {target} already exists with different content.")
                    overwrite = typer.confirm("Do you want to overwrite it?", default=False)
                    if overwrite:
                        shutil.copy2(tsv_path, target)
                        typer.echo(f"Overwrote {target}")
                    else:
                        raise typer.Abort()
            else:
                shutil.copy2(tsv_path, target)

            datasets.append(name)

        elif p.suffix.lower() == ".bed":

            if not p.exists():
                raise typer.BadParameter(f"BED file '{p}' does not exist.")

            bed_path = p.resolve()
            stem = bed_path.stem

            if stem.startswith("id_"):
                name = stem[3:]
            else:
                name = stem

            target = input_dir / f"id_{name}.bed"

            # Preprocess input file to a temp file first (filter + sort)
            with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".bed") as tmp_f:
                tmp_path = tmp_f.name
                # 1) Clean/filter rows
                with bed_path.open() as in_f:
                    for line in in_f:
                        line = line.rstrip("\n")
                        if not line:
                            continue  # skip empty lines

                        fields = line.split("\t")
                        if len(fields) < 4:
                            # Not enough columns to have chr, start, end, SVtype
                            continue

                        chrom = fields[0]
                        # Ensure chr prefix
                        if not chrom.startswith("chr"):
                            chrom = "chr" + chrom
                        fields[0] = chrom

                        # Keep only chr1-22, chrX, chrY
                        if chrom not in allowed_chroms:
                            continue

                        svtype = fields[3]
                        # Keep only DEL, DUP, INS, INV
                        if svtype not in allowed_svtypes:
                            continue

                        tmp_f.write("\t".join(fields) + "\n")

                tmp_f.flush()

            # 2) Sort cleaned rows into a second temp file
            with tempfile.NamedTemporaryFile(mode="w+", delete=False, suffix=".bed") as sorted_tmp:
                sorted_tmp_path = sorted_tmp.name
                subprocess.run(
                    ["sort", "-k1,1", "-k2,2n", tmp_path],
                    stdout=sorted_tmp,
                    check=True,
                )

            # Clean up first temp file
            os.unlink(tmp_path)

            # Now compare preprocessed content with existing target
            if target.exists():
                # Compare content of preprocessed file with existing target
                with open(sorted_tmp_path) as new_f:
                    new_content = new_f.read()
                with open(target) as existing_f:
                    existing_content = existing_f.read()

                if new_content == existing_content:
                    # Same content after preprocessing, just reuse existing
                    typer.echo(f"Using existing preprocessed file: {target}")
                    os.unlink(sorted_tmp_path)
                else:
                    # Different content, ask user what to do
                    typer.echo(f"File {target} already exists with different content.")
                    overwrite = typer.confirm("Do you want to overwrite it?", default=False)
                    if overwrite:
                        shutil.move(sorted_tmp_path, target)
                        typer.echo(f"Overwrote {target}")
                    else:
                        os.unlink(sorted_tmp_path)
                        raise typer.Abort()
            else:
                # Target doesn't exist, just move the preprocessed file
                shutil.move(sorted_tmp_path, target)

            datasets.append(name)

        else:

            name = raw
            if sequence_only:
                expected = input_dir / f"id_{name}.tsv"
            else:
                expected = input_dir / f"id_{name}.bed"
            if not expected.exists():
                raise typer.BadParameter(
                    f"'{name}' looks like a dataset name, but {expected} does not exist.\n"
                    f"Either:\n"
                    f"  - call CADD-SV with a BED/TSV file (recommended), e.g. 'caddsv {name}.bed'\n"
                    f"  - create {expected} manually."
                )
            datasets.append(name)

    datasets_value = ",".join(datasets)

    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(WORKFLOW_DIR / "Snakefile"),
        "--configfile",
        str(cfg),
        "--cores",
        str(threads),
        "--rerun-incomplete",
        "--use-conda",
        "--config",
        f"dataset={datasets_value}",
        f"mode={mode}",
        f"sequence_model={'True' if sequence_model else 'False'}",
    ]
    if force:
        cmd.append("--forceall")
    if unlock:
        cmd.append("--unlock")

    typer.echo("Running:\n  " + " ".join(cmd))

    env = os.environ.copy()
    env["PYTHONNOUSERSITE"] = "1"

    if check_time:
        start_wall = time.perf_counter()
        start_usage = resource.getrusage(resource.RUSAGE_CHILDREN)
        started_at = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        log_stamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        result = subprocess.run(cmd, env=env)

        end_wall = time.perf_counter()
        end_usage = resource.getrusage(resource.RUSAGE_CHILDREN)

        wall_time = end_wall - start_wall
        user_time = end_usage.ru_utime - start_usage.ru_utime
        sys_time = end_usage.ru_stime - start_usage.ru_stime
        cpu_time = user_time + sys_time
        cpu_pct = (cpu_time / wall_time * 100.0) if wall_time > 0 else 0.0

        log_lines = [
            f"[{started_at}] CADD-SV resource summary",
            f"Command: {' '.join(cmd)}",
            f"Return code: {result.returncode}",
            f"Wall time: {format_seconds(wall_time)}",
            f"CPU time: {format_seconds(cpu_time)} (user: {format_seconds(user_time)}, sys: {format_seconds(sys_time)})",
            f"CPU utilization: {cpu_pct:.1f}%",
            f"Max RSS (children, platform-dependent units): {end_usage.ru_maxrss}",
        ]

        for line in log_lines:
            typer.echo(line)

        log_path = Path(f"caddsv_run_{log_stamp}.log")
        log_path.write_text("\n".join(log_lines) + "\n")
        typer.echo(f"Resource log written to: {log_path.resolve()}")

        if result.returncode != 0:
            raise typer.Exit(code=result.returncode)
    else:
        subprocess.run(cmd, check=True, env=env)

    # Handle output based on mode
    if sequence_only:
        # Sequence-only mode: copy scored output to scored/ directory
        for name in datasets:
            src = Path("output") / f"{name}_seqonly_score.tsv"
            if not src.exists():
                typer.echo(f"Warning: Expected output {src} not found.")
                raise typer.Exit(code=1)

            dst = outdir / f"{name}_seqonly_score.tsv"

            if dst.exists():
                dst.unlink()

            shutil.copy2(src, dst)

        typer.echo(f"Sequence-only scores written to: {outdir.resolve()}")
    else:
        # Standard scoring mode
        for name in datasets:
            src = Path("output") / f"{name}bed_score100.bed"
            if not src.exists():
                raise typer.Exit(code=1)

            dst = outdir / f"{name}_score.tsv"

            if dst.exists():
                dst.unlink()

            # Always copy, no symlink (to be changed because redundant)
            shutil.copy2(src, dst)

        typer.echo(f"Final scores written to: {outdir.resolve()}")


def main():
    app()


if __name__ == "__main__":
    main()
