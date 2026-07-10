import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from typer.testing import CliRunner

from caddsv.cli import _build_snakemake_command, app


class BuildSnakemakeCommandTests(unittest.TestCase):
    def build_command(
        self,
        *,
        use_conda=True,
        conda_prefix=Path("/conda"),
        use_singularity=False,
        singularity_prefix=None,
        singularity_args=None,
    ):
        return _build_snakemake_command(
            cfg=Path("config.yml"),
            results_dir=Path("results"),
            threads=4,
            conda_prefix=conda_prefix,
            datasets_value="sample",
            mode="scoring",
            sequence_model=False,
            all_scores=False,
            annot_dir="/annotations",
            segmentnt_model_dir="/annotations/segment_nt",
            force=False,
            unlock=False,
            use_conda=use_conda,
            use_singularity=use_singularity,
            singularity_prefix=singularity_prefix,
            singularity_args=singularity_args,
        )

    def test_conda_is_enabled_by_default(self):
        command = self.build_command()

        self.assertEqual(command[:3], [sys.executable, "-m", "snakemake"])
        self.assertIn("--use-conda", command)
        prefix_index = command.index("--conda-prefix")
        self.assertEqual(command[prefix_index + 1], "/conda")
        frontend_index = command.index("--conda-frontend")
        self.assertEqual(command[frontend_index + 1], "conda")

    def test_conda_flags_are_omitted_when_disabled(self):
        command = self.build_command(use_conda=False, conda_prefix=None)

        self.assertNotIn("--use-conda", command)
        self.assertNotIn("--conda-prefix", command)
        self.assertIn("dataset=sample", command)

    def test_conda_prefix_is_required_when_conda_is_enabled(self):
        with self.assertRaisesRegex(ValueError, "conda_prefix is required"):
            self.build_command(conda_prefix=None)

    def test_singularity_backend_uses_image_cache_and_args(self):
        command = self.build_command(
            use_conda=False,
            conda_prefix=None,
            use_singularity=True,
            singularity_prefix=Path("/images"),
            singularity_args="--nv",
        )

        self.assertIn("--use-singularity", command)
        prefix_index = command.index("--singularity-prefix")
        self.assertEqual(command[prefix_index + 1], "/images")
        args_index = command.index("--singularity-args")
        self.assertEqual(command[args_index + 1], "--nv")
        self.assertNotIn("--use-conda", command)

    def test_conda_and_singularity_backends_cannot_be_combined(self):
        with self.assertRaisesRegex(ValueError, "mutually exclusive"):
            self.build_command(
                use_singularity=True,
                singularity_prefix=Path("/images"),
            )


class RunCommandTests(unittest.TestCase):
    def test_no_use_conda_uses_parent_environment_without_creating_cache(self):
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            results = root / "results"
            input_dir = results / "input"
            input_dir.mkdir(parents=True)
            (input_dir / "id_sample.bed").write_text("chr1\t1\t2\tDEL\n")

            workflow_output = results / "beds" / "sample" / "output"
            workflow_output.mkdir(parents=True)
            (workflow_output / "samplebed_score100.bed").write_text(
                "chr1\t1\t2\tDEL\n"
            )

            cache_home = root / "cache"
            with patch("caddsv.cli.subprocess.run") as run_mock:
                result = runner.invoke(
                    app,
                    [
                        "run",
                        "sample",
                        "--output-dir",
                        str(results),
                        "--no-use-conda",
                    ],
                    env={"XDG_CACHE_HOME": str(cache_home)},
                )

            self.assertEqual(result.exit_code, 0, result.output)
            command = run_mock.call_args.args[0]
            self.assertNotIn("--use-conda", command)
            self.assertNotIn("--conda-prefix", command)
            self.assertFalse(cache_home.exists())

    def test_singularity_mode_uses_its_cache_and_does_not_use_conda(self):
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            results = root / "results"
            input_dir = results / "input"
            input_dir.mkdir(parents=True)
            (input_dir / "id_sample.bed").write_text("chr1\t1\t2\tDEL\n")

            workflow_output = results / "beds" / "sample" / "output"
            workflow_output.mkdir(parents=True)
            (workflow_output / "samplebed_score100.bed").write_text(
                "chr1\t1\t2\tDEL\n"
            )

            cache_home = root / "cache"
            with patch("caddsv.cli.subprocess.run") as run_mock:
                result = runner.invoke(
                    app,
                    [
                        "run",
                        "sample",
                        "--output-dir",
                        str(results),
                        "--use-singularity",
                        "--singularity-args=--nv",
                    ],
                    env={"XDG_CACHE_HOME": str(cache_home)},
                )

            self.assertEqual(result.exit_code, 0, result.output)
            command = run_mock.call_args.args[0]
            self.assertIn("--use-singularity", command)
            self.assertNotIn("--use-conda", command)
            self.assertIn("--singularity-args", command)
            args_index = command.index("--singularity-args")
            self.assertIn("--nv", command[args_index + 1])
            self.assertTrue((cache_home / "caddsv" / "snakemake-singularity").exists())


if __name__ == "__main__":
    unittest.main()
