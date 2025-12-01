from pathlib import Path
import subprocess
import sys
from typing import List, Optional
import typer
import os
import shutil

app = typer.Typer(add_completion=False, help="CADD-SV Snakemake-based scoring tool")

PKG_DIR = Path(__file__).resolve().parent
WORKFLOW_DIR = PKG_DIR / "workflow"
DEFAULT_CONFIG = PKG_DIR / "config.yml"


def _cli(
    items: List[str] = typer.Argument(
        ...,
        help=(
            "One or more inputs, each either:\n"
            "  - a BED file (e.g. variants.bed), or\n"
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
):
    """
    Single-command interface: call as

        caddsv my_calls.bed
        caddsv input/id_longrange_cherie.bed
        caddsv longrange_cherie
    """
    cfg = config if config is not None else DEFAULT_CONFIG

    # Where Snakemake expects id_<name>.bed
    input_dir = Path("input")
    input_dir.mkdir(exist_ok=True)

    datasets: List[str] = []

    for raw in items:
        p = Path(raw)

        if p.suffix.lower() == ".bed":
            # Case 1: user gave a BED path (absolute or relative)
            if not p.exists():
                raise typer.BadParameter(f"BED file '{p}' does not exist.")

            bed_path = p.resolve()
            stem = bed_path.stem

            # If the file is already named id_<name>.bed, strip the id_ prefix
            if stem.startswith("id_"):
                name = stem[3:]
            else:
                name = stem

            target = input_dir / f"id_{name}.bed"

            # If target exists and points elsewhere, complain
            if target.exists() and target.resolve() != bed_path:
                raise typer.BadParameter(
                    f"Input conflict: {target} already exists and does not match {bed_path}."
                )

            # Create symlink or copy if needed
            if not target.exists():
                try:
                    target.symlink_to(bed_path)
                except OSError:
                    shutil.copy2(bed_path, target)

            datasets.append(name)

        else:
            # Case 2: user gave a dataset name
            name = raw
            expected = input_dir / f"id_{name}.bed"
            if not expected.exists():
                raise typer.BadParameter(
                    f"'{name}' looks like a dataset name, but {expected} does not exist.\n"
                    f"Either:\n"
                    f"  - call CADD-SV with a BED file (recommended), e.g. 'caddsv {name}.bed'\n"
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
    ]
    if force:
        cmd.append("--forceall")

    typer.echo("Running:\n  " + " ".join(cmd))

    env = os.environ.copy()
    env["PYTHONNOUSERSITE"] = "1"

    subprocess.run(cmd, check=True, env=env)


def main():
    # Entry point used by the console script:
    # caddsv -> this function -> typer.run(_cli)
    typer.run(_cli)


if __name__ == "__main__":
    main()
