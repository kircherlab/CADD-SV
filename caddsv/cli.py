from pathlib import Path, PurePosixPath
from decimal import Decimal, InvalidOperation, ROUND_HALF_UP
import hashlib
from importlib.metadata import PackageNotFoundError, version as package_version
import shlex
import subprocess
import sys
from typing import List, Optional
import tempfile
import tarfile
import typer
import os
import shutil
import time
from datetime import datetime
import resource
import urllib.error
import urllib.request

from caddsv.container_images import (
    COORDINATE_BASED_CONTAINER_ENVIRONMENTS,
    DEFAULT_CONTAINER_IMAGES,
)

app = typer.Typer(add_completion=False, help="CADD-SV Snakemake-based scoring tool")


def _version_callback(value: bool) -> None:
    if not value:
        return

    try:
        installed_version = package_version("caddsv")
    except PackageNotFoundError:
        installed_version = "unknown (package metadata unavailable)"

    typer.echo(installed_version)
    raise typer.Exit()


@app.callback()
def cli_callback(
    version: Optional[bool] = typer.Option(
        None,
        "--version",
        callback=_version_callback,
        is_eager=True,
        help="Show the installed CADD-SV version and exit.",
    ),
) -> None:
    """CADD-SV command-line interface."""

PKG_DIR = Path(__file__).resolve().parent
WORKFLOW_DIR = PKG_DIR / "workflow"
MODELS_DIR = WORKFLOW_DIR / "models"
DEFAULT_CONFIG = PKG_DIR / "config.yml"
CONDA_ENVS_SNAKEFILE = WORKFLOW_DIR / "create_envs.smk"
CONDA_PREFIX_ENV_VAR = "CADD_SV_CONDA_PREFIX"
CONDA_CACHE_SUBDIR = Path("caddsv") / "snakemake-conda"
SINGULARITY_PREFIX_ENV_VARS = (
    "CADD_SV_SINGULARITY_PREFIX",
    "APPTAINER_CACHEDIR",
    "SINGULARITY_CACHEDIR",
)
SINGULARITY_CACHE_SUBDIR = Path("caddsv") / "snakemake-singularity"
SEGMENTNT_REPO_ID = "InstaDeepAI/segment_nt"
SEGMENTNT_DIRNAME = "segment_nt"
ANNOTATIONS_ARCHIVE_URL = "https://kircherlab.bihealth.org/download/CADD-SV/v2.0/dependencies.tar.gz"
DOWNLOAD_CHUNK_SIZE = 8 * 1024 * 1024
SEGMENTNT_ALLOW_PATTERNS = [
    "README.md",
    "config.json",
    "modeling_segment_nt.py",
    "pytorch_model.bin",
    "segment_nt_config.py",
    "special_tokens_map.json",
    "tokenizer_config.json",
    "vocab.txt",
]
SEGMENTNT_NOTICE = """SegmentNT model files
Source: https://huggingface.co/InstaDeepAI/segment_nt
Developed by: InstaDeep
License: Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)
License URL: https://creativecommons.org/licenses/by-nc-sa/4.0/

These files are downloaded from Hugging Face for local CADD-SV SegmentNT
inference. They are not CADD-SV-authored files. Commercial use and
redistribution are subject to the upstream SegmentNT license.
"""

_RAW_SCORE_SUFFIX = "_score"
_SCORE_PRECISION = Decimal("0.0001")


def _format_raw_score(value: str) -> str:
    """Render a finite raw score with exactly four decimal places."""
    try:
        score = Decimal(value)
    except (InvalidOperation, ValueError):
        return value

    if not score.is_finite():
        return value

    return format(score.quantize(_SCORE_PRECISION, rounding=ROUND_HALF_UP), "f")


def write_final_score_file(source: Path, destination: Path) -> None:
    """Create the user-facing score file without changing workflow output."""
    with source.open(encoding="utf-8") as input_file:
        header_line = input_file.readline()
        if not header_line:
            destination.write_text("", encoding="utf-8")
            return

        header = header_line.rstrip("\r\n").split("\t")
        if header[0] == "chr":
            header[0] = "#chr"
        raw_score_indices = [
            index
            for index, column in enumerate(header)
            if column.endswith(_RAW_SCORE_SUFFIX)
        ]

        with destination.open("w", encoding="utf-8", newline="\n") as output_file:
            output_file.write("\t".join(header) + "\n")
            for line in input_file:
                fields = line.rstrip("\r\n").split("\t")
                for index in raw_score_indices:
                    if index < len(fields):
                        fields[index] = _format_raw_score(fields[index])
                output_file.write("\t".join(fields) + "\n")


def _download_file(url: str, destination: Path) -> None:
    try:
        with urllib.request.urlopen(url, timeout=60) as response:
            with destination.open("wb") as output:
                shutil.copyfileobj(response, output, length=DOWNLOAD_CHUNK_SIZE)
    except urllib.error.URLError as exc:
        raise typer.BadParameter(f"Failed to download {url}: {exc}") from exc


def _strip_archive_root(member_name: str) -> Optional[PurePosixPath]:
    member_path = PurePosixPath(member_name)
    if member_path.is_absolute():
        raise typer.BadParameter(
            f"Refusing to extract absolute archive path: {member_name}"
        )

    parts = [part for part in member_path.parts if part not in ("", ".")]
    if any(part == ".." for part in parts):
        raise typer.BadParameter(
            f"Refusing to extract unsafe archive path: {member_name}"
        )

    stripped_parts = parts[1:]
    if not stripped_parts:
        return None
    return PurePosixPath(*stripped_parts)


def _ensure_within_directory(base: Path, destination: Path) -> None:
    base = base.resolve()
    destination = destination.resolve()
    try:
        destination.relative_to(base)
    except ValueError as exc:
        raise typer.BadParameter(
            f"Refusing to extract outside annotation directory: {destination}"
        ) from exc


def _extract_annotations_archive(archive_path: Path, target: Path) -> None:
    target.mkdir(parents=True, exist_ok=True)
    with tarfile.open(archive_path, "r:*") as archive:
        for member in archive:
            relative_path = _strip_archive_root(member.name)
            if relative_path is None:
                continue

            destination = target / Path(*relative_path.parts)
            _ensure_within_directory(target, destination)

            if member.isdir():
                destination.mkdir(parents=True, exist_ok=True)
            elif member.isfile():
                destination.parent.mkdir(parents=True, exist_ok=True)
                source = archive.extractfile(member)
                if source is None:
                    raise typer.BadParameter(
                        f"Could not read archive entry: {member.name}"
                    )
                with source, destination.open("wb") as output:
                    shutil.copyfileobj(source, output, length=DOWNLOAD_CHUNK_SIZE)
                os.chmod(destination, member.mode & 0o777)
            else:
                raise typer.BadParameter(
                    f"Refusing to extract unsupported archive entry: {member.name}"
                )


def download_annotations(target: Path, url: str = ANNOTATIONS_ARCHIVE_URL) -> None:
    target.mkdir(parents=True, exist_ok=True)
    archive_path: Optional[Path] = None

    with tempfile.NamedTemporaryFile(
        prefix="caddsv-dependencies-",
        suffix=".tar.gz",
        dir=target,
        delete=False,
    ) as archive_file:
        archive_path = Path(archive_file.name)

    try:
        typer.echo(f"Downloading dependencies from: {url}")
        _download_file(url, archive_path)
        typer.echo("Uncompressing dependencies...")
        _extract_annotations_archive(archive_path, target)
    finally:
        if archive_path is not None:
            archive_path.unlink(missing_ok=True)


def write_segmentnt_notice(target: Path) -> None:
    (target / "SEGMENTNT_NOTICE.txt").write_text(SEGMENTNT_NOTICE)


def download_segmentnt(target: Path, repo_id: str, force: bool = False) -> None:
    try:
        from huggingface_hub import snapshot_download
    except ImportError as exc:
        raise typer.BadParameter(
            "Downloading SegmentNT requires huggingface_hub. "
            "Install CADD-SV with its current dependencies or run "
            "'pip install huggingface_hub'."
        ) from exc

    target.mkdir(parents=True, exist_ok=True)
    if force:
        for path in target.iterdir():
            if path.is_dir():
                shutil.rmtree(path)
            else:
                path.unlink()

    typer.echo(f"Downloading SegmentNT from {repo_id} to: {target}")
    snapshot_download(
        repo_id=repo_id,
        local_dir=str(target),
        allow_patterns=SEGMENTNT_ALLOW_PATTERNS,
    )
    write_segmentnt_notice(target)
    typer.echo(f"SegmentNT model files are available at: {target}")


def format_seconds(value: float) -> str:
    return f"{value:.2f}s"


def _default_conda_prefix() -> Path:
    cache_home = os.environ.get("XDG_CACHE_HOME")
    cache_base = (
        Path(cache_home).expanduser()
        if cache_home
        else Path.home() / ".cache"
    )
    return cache_base / CONDA_CACHE_SUBDIR


def resolve_conda_prefix(conda_prefix: Optional[Path]) -> Path:
    if conda_prefix is not None:
        return Path(conda_prefix).expanduser().resolve()

    env_value = os.environ.get(CONDA_PREFIX_ENV_VAR)
    if env_value:
        return Path(env_value).expanduser().resolve()

    return _default_conda_prefix().expanduser().resolve()


def _default_singularity_prefix() -> Path:
    cache_home = os.environ.get("XDG_CACHE_HOME")
    cache_base = (
        Path(cache_home).expanduser()
        if cache_home
        else Path.home() / ".cache"
    )
    return cache_base / SINGULARITY_CACHE_SUBDIR


def resolve_singularity_prefix(singularity_prefix: Optional[Path]) -> Path:
    if singularity_prefix is not None:
        return Path(singularity_prefix).expanduser().resolve()

    for env_var in SINGULARITY_PREFIX_ENV_VARS:
        env_value = os.environ.get(env_var)
        if env_value:
            return Path(env_value).expanduser().resolve()

    return _default_singularity_prefix().expanduser().resolve()


def container_image_path(image: str, prefix: Path) -> Path:
    """Return the cache path Snakemake derives from a container URI."""
    image_hash = hashlib.md5(image.encode(), usedforsecurity=False).hexdigest()
    return prefix / f"{image_hash}.simg"


def find_container_runtime() -> str:
    for executable in ("apptainer", "singularity"):
        runtime = shutil.which(executable)
        if runtime:
            return runtime
    raise typer.BadParameter(
        "Downloading environment images requires Apptainer or Singularity. "
        "Install one of them and ensure its executable is available on PATH."
    )


def download_container_images(
    target: Path,
    coordinate_based_only: bool = False,
    force: bool = False,
) -> None:
    runtime = find_container_runtime()
    target.mkdir(parents=True, exist_ok=True)
    environments = (
        COORDINATE_BASED_CONTAINER_ENVIRONMENTS
        if coordinate_based_only
        else tuple(DEFAULT_CONTAINER_IMAGES)
    )

    for environment in environments:
        image = DEFAULT_CONTAINER_IMAGES[environment]
        destination = container_image_path(image, target)

        if destination.exists() and not force:
            typer.echo(f"Using cached {environment} image: {destination}")
            continue

        action = "Refreshing" if destination.exists() else "Downloading"
        typer.echo(f"{action} {environment} image from: {image}")
        temporary_dir = Path(
            tempfile.mkdtemp(prefix=".caddsv-image-", dir=target)
        )
        temporary_image = temporary_dir / destination.name
        try:
            subprocess.run(
                [runtime, "pull", str(temporary_image), image],
                check=True,
            )
            if not temporary_image.is_file():
                raise typer.BadParameter(
                    f"The container runtime did not create the expected image: "
                    f"{temporary_image}"
                )
            os.replace(temporary_image, destination)
        except subprocess.CalledProcessError as exc:
            raise typer.BadParameter(
                f"Failed to download {environment} image from {image}."
            ) from exc
        finally:
            shutil.rmtree(temporary_dir, ignore_errors=True)

        typer.echo(f"Cached {environment} image at: {destination}")

    typer.echo(f"Environment images are available at: {target}")


def build_conda_envs_command(
    target: Path,
    workdir: Path,
    coordinate_based_only: bool = False,
) -> List[str]:
    workflow_target = (
        "coordinate_based_envs" if coordinate_based_only else "all_envs"
    )
    return [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(CONDA_ENVS_SNAKEFILE),
        "--directory",
        str(workdir),
        "--cores",
        "1",
        "--use-conda",
        "--conda-prefix",
        str(target),
        "--conda-frontend",
        "conda",
        "--conda-create-envs-only",
        workflow_target,
    ]


def create_conda_environments(
    target: Path,
    coordinate_based_only: bool = False,
) -> None:
    target.mkdir(parents=True, exist_ok=True)
    typer.echo(f"Creating Conda environments in: {target}")

    with tempfile.TemporaryDirectory(prefix="caddsv-conda-envs-") as workdir:
        command = build_conda_envs_command(
            target,
            Path(workdir),
            coordinate_based_only=coordinate_based_only,
        )
        try:
            subprocess.run(command, check=True)
        except subprocess.CalledProcessError as exc:
            raise typer.BadParameter(
                f"Failed to create Conda environments in {target}."
            ) from exc

    typer.echo(f"Conda environments are available at: {target}")


def build_singularity_args(
    annot_dir: str,
    segmentnt_model_dir: str,
    user_args: Optional[str],
) -> str:
    """Bind packaged scripts and external data at their absolute host paths."""
    bind_paths = [WORKFLOW_DIR, Path(annot_dir), Path(segmentnt_model_dir)]
    seen = set()
    args: List[str] = []

    for path in bind_paths:
        resolved = path.expanduser().resolve()
        if not resolved.exists() or resolved in seen:
            continue
        seen.add(resolved)
        binding = f"{resolved}:{resolved}:ro"
        args.extend(["--bind", shlex.quote(binding)])

    if user_args:
        args.append(user_args)

    return " ".join(args)


def parse_snakemake_args(value: Optional[str]) -> List[str]:
    """Parse optional user-provided Snakemake arguments safely."""
    if not value:
        return []

    try:
        return shlex.split(value)
    except ValueError as exc:
        raise typer.BadParameter(
            f"Invalid --snakemake-args value: {exc}"
        ) from exc


def resolve_scaler_dir(scaler_dir: Optional[Path]) -> Path:
    if scaler_dir is None:
        raise typer.BadParameter(
            "caddsv scale requires --scaler-dir pointing to a directory "
            "containing {TYPE}_stats.tsv files."
        )

    resolved = Path(scaler_dir).expanduser().resolve()
    if not resolved.exists():
        raise typer.BadParameter(f"Scaler directory does not exist: {resolved}")
    if not resolved.is_dir():
        raise typer.BadParameter(f"Scaler path is not a directory: {resolved}")
    return resolved


def _build_snakemake_command(
    cfg: Path,
    results_dir: Path,
    threads: int,
    conda_prefix: Optional[Path],
    datasets_value: str,
    mode: str,
    sequence_model: bool,
    all_scores: bool,
    annot_dir: str,
    segmentnt_model_dir: str,
    force: bool,
    unlock: bool,
    use_conda: bool = True,
    use_singularity: bool = False,
    singularity_prefix: Optional[Path] = None,
    singularity_args: Optional[str] = None,
    extra_snakemake_args: Optional[List[str]] = None,
) -> List[str]:
    if use_conda and use_singularity:
        raise ValueError("conda and singularity backends are mutually exclusive")

    cmd = [
        sys.executable,
        "-m",
        "snakemake",
        "--snakefile",
        str(WORKFLOW_DIR / "Snakefile"),
        "--configfile",
        str(cfg),
        "--directory",
        str(results_dir),
        "--cores",
        str(threads),
        "--rerun-incomplete",
    ]
    if use_conda:
        if conda_prefix is None:
            raise ValueError("conda_prefix is required when use_conda is enabled")
        cmd.extend(
            [
                "--use-conda",
                "--conda-prefix",
                str(conda_prefix),
                "--conda-frontend",
                "conda",
            ]
        )
    elif use_singularity:
        if singularity_prefix is None:
            raise ValueError(
                "singularity_prefix is required when singularity is enabled"
            )
        cmd.extend(
            [
                "--use-singularity",
                "--singularity-prefix",
                str(singularity_prefix),
            ]
        )
        if singularity_args:
            cmd.extend(["--singularity-args", singularity_args])
    cmd.extend(
        [
            "--config",
            f"dataset={datasets_value}",
            f"mode={mode}",
            f"sequence_model={'True' if sequence_model else 'False'}",
            f"all_scores={'True' if all_scores else 'False'}",
            f"annotations_dir={annot_dir}",
            f"segmentnt_model_dir={segmentnt_model_dir}",
        ]
    )
    if force:
        cmd.append("--forceall")
    if unlock:
        cmd.append("--unlock")
    if extra_snakemake_args:
        cmd.extend(extra_snakemake_args)
    return cmd


@app.command()
def get(
        flag: str = typer.Argument(None),
        annotations_dir: Optional[Path] = typer.Option(
            None, "--annotations-dir",
            help="Directory to store annotations (default: ./annotations)"
        ),
        with_segmentnt: bool = typer.Option(
            False,
            "--with-segmentnt",
            help="Also download SegmentNT model files into <annotations-dir>/segment_nt.",
        ),
        force_segmentnt: bool = typer.Option(
            False,
            "--force-segmentnt",
            help="Replace an existing local SegmentNT model directory.",
        ),
        segmentnt_repo: str = typer.Option(
            SEGMENTNT_REPO_ID,
            "--segmentnt-repo",
            help="Hugging Face repository to use for SegmentNT model files.",
        ),
        singularity_prefix: Optional[Path] = typer.Option(
            None,
            "--singularity-prefix",
            "--apptainer-prefix",
            help=(
                "Directory to store prefetched environment images "
                "(default: the CADD-SV Singularity/Apptainer cache)."
            ),
        ),
        conda_prefix: Optional[Path] = typer.Option(
            None,
            "--conda-prefix",
            help=(
                "Directory for prefetched Snakemake Conda environments "
                "(default: CADD_SV_CONDA_PREFIX or user cache)."
            ),
        ),
        use_conda: bool = typer.Option(
            False,
            "--use-conda",
            help="Create Conda environments instead of downloading images.",
        ),
        coordinate_based_only: bool = typer.Option(
            False,
            "--coordinate-based-only",
            help="Prepare preprocessing, SV, and training environments, but not NT.",
        ),
        force_envs: bool = typer.Option(
            False,
            "--force-envs",
            help="Replace environment images already present in the image cache.",
        ),
):

    target = Path(annotations_dir).resolve() if annotations_dir else Path("annotations").resolve()

    if flag == "annotations":
        download_annotations(target)
        typer.echo(f"Annotations extracted to: {target}")
        if with_segmentnt:
            download_segmentnt(target / SEGMENTNT_DIRNAME, segmentnt_repo, force_segmentnt)
        typer.echo("DONE")
    elif flag in {"segmentnt", "segmentNT", "SegmentNT"}:
        download_segmentnt(target / SEGMENTNT_DIRNAME, segmentnt_repo, force_segmentnt)
        typer.echo("DONE")
    elif flag == "envs":
        if use_conda:
            if singularity_prefix is not None:
                raise typer.BadParameter(
                    "--singularity-prefix/--apptainer-prefix cannot be used "
                    "with --use-conda."
                )
            if force_envs:
                raise typer.BadParameter(
                    "--force-envs only applies to environment images. Conda "
                    "environments are recreated automatically when their YAML "
                    "definitions change."
                )
            conda_target = resolve_conda_prefix(conda_prefix)
            create_conda_environments(
                conda_target,
                coordinate_based_only=coordinate_based_only,
            )
        else:
            if conda_prefix is not None:
                raise typer.BadParameter(
                    "--conda-prefix requires --use-conda."
                )
            image_target = resolve_singularity_prefix(singularity_prefix)
            download_container_images(
                image_target,
                coordinate_based_only=coordinate_based_only,
                force=force_envs,
            )
        typer.echo("DONE")
    else:
        raise typer.BadParameter("Expected 'annotations', 'segmentnt', or 'envs'.")

@app.command()
def run(
    ctx: typer.Context,
    items: List[str] = typer.Argument(
        ...,
        help=(
            "One or more inputs, each either:\n"
            "  - a BED file (e.g. variants.bed), or\n"
            "  - a TSV file with --seqonly (REF, ALT, [TYPE], [ID]), or\n"
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
    sequence_model: bool = typer.Option(False, "--seqresolved", help="Optional to run the model with the integration of SegmentNT derived annotations"),
    sequence_only: bool = typer.Option(
        False, "--seqonly",
        help="Sequence-only mode: input TSV with columns REF, ALT, [TYPE], [ID]. "
             "Runs only SegmentNT, generates SB/SBref/DB features without coordinate-based annotations."
    ),
    all_scores: bool = typer.Option(
        False, "--allscores",
        hidden=True,
        help="Run all scoring models: default, sequence-resolved, and sequence-only."
    ),
    unlock: bool = typer.Option(
        False, "--unlock", help="snakemake --unlock"
    ),
    annotations_dir: Optional[Path] = typer.Option(
        None, "--annotations-dir",
        help="Path to annotation directory (default: ./annotations)"
    ),
    output_dir: Optional[Path] = typer.Option(
        None, "--output-dir", "-o",
        help="Results directory (default: ./caddsv_results)"
    ),
    conda_prefix: Optional[Path] = typer.Option(
        None,
        "--conda-prefix",
        help=(
            "Directory for Snakemake conda environments "
            "(default: CADD_SV_CONDA_PREFIX or user cache)."
        ),
    ),
    use_conda: bool = typer.Option(
        True,
        "--use-conda/--no-use-conda",
        help=(
            "Use Snakemake conda environments. Disable this when all workflow "
            "dependencies are already available in the parent environment."
        ),
    ),
    use_singularity: bool = typer.Option(
        False,
        "--use-singularity",
        "--use-apptainer",
        help=(
            "Run rules in their versioned Singularity/Apptainer containers "
            "instead of creating conda environments."
        ),
    ),
    singularity_prefix: Optional[Path] = typer.Option(
        None,
        "--singularity-prefix",
        "--apptainer-prefix",
        help=(
            "Directory for cached Singularity/Apptainer images "
            "(default: CADD_SV_SINGULARITY_PREFIX, runtime cache, or user cache)."
        ),
    ),
    singularity_args: Optional[str] = typer.Option(
        None,
        "--singularity-args",
        "--apptainer-args",
        help=(
            "Additional arguments passed to Singularity/Apptainer; use '--nv' "
            "for NVIDIA GPU access."
        ),
    ),
    check_time: bool = typer.Option(
        False,
        "--check-time",
        help="Track time and resource usage; log to a .log file.",
    ),
    snakemake_args: Optional[str] = typer.Option(
        None,
        "--snakemake-args",
        help=(
            "Additional Snakemake arguments. Parsed with shell-style quoting "
            "and appended last, so they can override CADD-SV options."
        ),
    ),
):
    cfg = config if config is not None else DEFAULT_CONFIG
    annot_dir = str(Path(annotations_dir).resolve()) if annotations_dir else str(Path("annotations").resolve())

    use_conda_was_explicit = (
        ctx.get_parameter_source("use_conda").name == "COMMANDLINE"
    )
    if use_singularity and use_conda and use_conda_was_explicit:
        raise typer.BadParameter(
            "--use-conda and --use-singularity are mutually exclusive."
        )
    if use_singularity:
        use_conda = False

    if not use_conda and conda_prefix is not None:
        raise typer.BadParameter(
            "--conda-prefix can only be used with the conda backend."
        )
    if not use_singularity and singularity_prefix is not None:
        raise typer.BadParameter(
            "--singularity-prefix requires --use-singularity."
        )
    if not use_singularity and singularity_args is not None:
        raise typer.BadParameter(
            "--singularity-args requires --use-singularity."
        )

    extra_snakemake_args = parse_snakemake_args(snakemake_args)

    # Override mode if sequence_only flag is set
    if sequence_only:
        mode = "seqonly"

    # all_scores implies sequence_model
    if all_scores:
        sequence_model = True

    results_dir = Path(output_dir) if output_dir else Path("caddsv_results")
    results_dir.mkdir(exist_ok=True)

    outdir = results_dir / "scored"
    outdir.mkdir(exist_ok=True)

    input_dir = results_dir / "input"
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
    resolved_conda_prefix: Optional[Path] = None
    resolved_singularity_prefix: Optional[Path] = None
    resolved_singularity_args: Optional[str] = None
    if use_conda:
        resolved_conda_prefix = resolve_conda_prefix(conda_prefix)
        resolved_conda_prefix.mkdir(parents=True, exist_ok=True)
    segmentnt_model_dir = os.environ.get(
        "SEGMENTNT_MODEL",
        str(Path(annot_dir) / SEGMENTNT_DIRNAME),
    )
    if use_singularity:
        resolved_singularity_prefix = resolve_singularity_prefix(
            singularity_prefix
        )
        resolved_singularity_prefix.mkdir(parents=True, exist_ok=True)
        resolved_singularity_args = build_singularity_args(
            annot_dir,
            segmentnt_model_dir,
            singularity_args,
        )
    cmd = _build_snakemake_command(
        cfg=Path(cfg),
        results_dir=results_dir,
        threads=threads,
        conda_prefix=resolved_conda_prefix,
        datasets_value=datasets_value,
        mode=mode,
        sequence_model=sequence_model,
        all_scores=all_scores,
        annot_dir=annot_dir,
        segmentnt_model_dir=segmentnt_model_dir,
        force=force,
        unlock=unlock,
        use_conda=use_conda,
        use_singularity=use_singularity,
        singularity_prefix=resolved_singularity_prefix,
        singularity_args=resolved_singularity_args,
        extra_snakemake_args=extra_snakemake_args,
    )

    typer.echo("Running:\n  " + " ".join(cmd))

    env = os.environ.copy()

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
            src = results_dir / "beds" / name / "output" / f"{name}_seqonly_score.tsv"
            if not src.exists():
                typer.echo(f"Warning: Expected output {src} not found.")
                raise typer.Exit(code=1)

            dst = outdir / f"{name}_seqonly_score.tsv"

            write_final_score_file(src, dst)

        typer.echo(f"Sequence-only scores written to: {outdir.resolve()}")
    else:
        # Standard scoring mode
        for name in datasets:
            src = results_dir / "beds" / name / "output" / f"{name}bed_score100.bed"
            if not src.exists():
                raise typer.Exit(code=1)

            dst = outdir / f"{name}_score.tsv"

            write_final_score_file(src, dst)

        typer.echo(f"Final scores written to: {outdir.resolve()}")


@app.command(hidden=True)
def scale(
    items: List[str] = typer.Argument(
        ...,
        help="One or more *_score.tsv files to z-score scale",
    ),
    output_dir: Optional[Path] = typer.Option(
        None, "--output-dir", "-o",
        help="Directory for scaled output (default: same directory as input)",
    ),
    scaler_dir: Optional[Path] = typer.Option(
        None, "--scaler-dir", "--stats-dir",
        help=(
            "Directory containing {TYPE}_stats.tsv scaler files. "
            "--stats-dir is accepted as a compatibility alias."
        ),
    ),
):
    """Generate z-score scaled features from scored CADD-SV output files."""
    sdir = resolve_scaler_dir(scaler_dir)
    from caddsv.workflow.scripts.scale_features import scale_features

    for raw in items:
        scored = Path(raw)
        if not scored.exists():
            raise typer.BadParameter(f"File not found: {scored}")

        stem = scored.stem
        # Replace _score suffix with _scaled, otherwise just append _scaled
        if stem.endswith("_score"):
            out_stem = stem[: -len("_score")] + "_scaled"
        else:
            out_stem = stem + "_scaled"

        if output_dir:
            dest = Path(output_dir)
            dest.mkdir(parents=True, exist_ok=True)
        else:
            dest = scored.parent

        out_path = dest / f"{out_stem}.tsv"
        scale_features(str(scored), str(out_path), sdir)
        typer.echo(f"Scaled features written to: {out_path}")


def main():
    app()


if __name__ == "__main__":
    main()
