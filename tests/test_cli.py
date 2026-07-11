import subprocess
import sys
import tempfile
import unittest
from importlib.metadata import version as package_version
from pathlib import Path
from unittest.mock import patch

from typer.testing import CliRunner

from caddsv.cli import (
    DEFAULT_CONFIG,
    _build_snakemake_command,
    app,
    container_image_path,
    write_final_score_file,
)
from caddsv.container_images import (
    COORDINATE_BASED_CONTAINER_ENVIRONMENTS,
    DEFAULT_CONTAINER_IMAGES,
)


class VersionTests(unittest.TestCase):
    def test_version_reports_installed_distribution_metadata(self):
        result = CliRunner().invoke(app, ["--version"])

        self.assertEqual(result.exit_code, 0, result.output)
        self.assertEqual(result.output.strip(), package_version("caddsv"))


class FinalScoreFileTests(unittest.TestCase):
    def test_final_file_formats_raw_scores_and_preserves_other_values(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "workflow-score.bed"
            destination = Path(tmpdir) / "final-score.tsv"
            source_contents = (
                "chr\tstart\tend\ttype\tCADD-SV_PHRED\tCADD-SV_score"
                "\tCADD-SV-SR_score\tfeature\n"
                "chr1\t1\t2\tDEL\t12.345678\t0.123456\t0.50005\t9.876543\n"
                "chr2\t3\t4\tDUP\tNaN\tNaN\tnot-a-number\t\n"
            )
            source.write_text(source_contents)

            write_final_score_file(source, destination)

            self.assertEqual(source.read_text(), source_contents)
            self.assertEqual(
                destination.read_text(),
                (
                    "#chr\tstart\tend\ttype\tCADD-SV_PHRED\tCADD-SV_score"
                    "\tCADD-SV-SR_score\tfeature\n"
                    "chr1\t1\t2\tDEL\t12.345678\t0.1235\t0.5001\t9.876543\n"
                    "chr2\t3\t4\tDUP\tNaN\tNaN\tnot-a-number\t\n"
                ),
            )

    def test_sequence_only_final_file_keeps_its_non_chromosome_header(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "workflow-score.tsv"
            destination = Path(tmpdir) / "final-score.tsv"
            source.write_text(
                "id\ttype\tCADD-SV_seqonly_PHRED\tCADD-SV_seqonly_score\n"
                "sample-1\tDEL\t7.123456\t0.5\n"
            )

            write_final_score_file(source, destination)

            self.assertEqual(
                destination.read_text(),
                "id\ttype\tCADD-SV_seqonly_PHRED\tCADD-SV_seqonly_score\n"
                "sample-1\tDEL\t7.123456\t0.5000\n",
            )


class GetEnvironmentImagesTests(unittest.TestCase):
    def setUp(self):
        self.runner = CliRunner()

    def write_fake_image(self, command, check):
        self.assertTrue(check)
        self.assertEqual(command[1], "pull")
        Path(command[2]).write_bytes(b"SIF")

    def test_packaged_config_matches_shared_image_defaults(self):
        config_text = DEFAULT_CONFIG.read_text()
        for environment, image in DEFAULT_CONTAINER_IMAGES.items():
            self.assertIn(f"{environment}: {image}", config_text)

    def test_get_envs_downloads_all_images_with_snakemake_names(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "images"
            with (
                patch("caddsv.cli.shutil.which", return_value="/usr/bin/apptainer"),
                patch(
                    "caddsv.cli.subprocess.run",
                    side_effect=self.write_fake_image,
                ) as run_mock,
            ):
                result = self.runner.invoke(
                    app,
                    ["get", "envs", "--apptainer-prefix", str(prefix)],
                )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(run_mock.call_count, 4)
            for call, image in zip(
                run_mock.call_args_list,
                DEFAULT_CONTAINER_IMAGES.values(),
            ):
                command = call.args[0]
                self.assertEqual(command[0], "/usr/bin/apptainer")
                self.assertEqual(command[3], image)
                self.assertTrue(container_image_path(image, prefix).is_file())
            self.assertFalse(list(prefix.glob(".caddsv-image-*")))

    def test_coordinate_based_only_omits_nt_image(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "images"
            with (
                patch("caddsv.cli.shutil.which", return_value="/usr/bin/apptainer"),
                patch(
                    "caddsv.cli.subprocess.run",
                    side_effect=self.write_fake_image,
                ) as run_mock,
            ):
                result = self.runner.invoke(
                    app,
                    [
                        "get",
                        "envs",
                        "--coordinate-based-only",
                        "--apptainer-prefix",
                        str(prefix),
                    ],
                )

            self.assertEqual(result.exit_code, 0, result.output)
            pulled_images = [call.args[0][3] for call in run_mock.call_args_list]
            self.assertEqual(
                pulled_images,
                [
                    DEFAULT_CONTAINER_IMAGES[environment]
                    for environment in COORDINATE_BASED_CONTAINER_ENVIRONMENTS
                ],
            )
            self.assertNotIn(DEFAULT_CONTAINER_IMAGES["nt"], pulled_images)

    def test_existing_images_are_skipped_and_force_envs_refreshes_them(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "images"
            prefix.mkdir()
            for image in DEFAULT_CONTAINER_IMAGES.values():
                container_image_path(image, prefix).write_bytes(b"old")

            with (
                patch("caddsv.cli.shutil.which", return_value="/usr/bin/apptainer"),
                patch("caddsv.cli.subprocess.run") as run_mock,
            ):
                result = self.runner.invoke(
                    app,
                    ["get", "envs", "--apptainer-prefix", str(prefix)],
                )

            self.assertEqual(result.exit_code, 0, result.output)
            run_mock.assert_not_called()

            with (
                patch("caddsv.cli.shutil.which", return_value="/usr/bin/apptainer"),
                patch(
                    "caddsv.cli.subprocess.run",
                    side_effect=self.write_fake_image,
                ) as run_mock,
            ):
                result = self.runner.invoke(
                    app,
                    [
                        "get",
                        "envs",
                        "--apptainer-prefix",
                        str(prefix),
                        "--force-envs",
                    ],
                )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(run_mock.call_count, 4)
            for image in DEFAULT_CONTAINER_IMAGES.values():
                self.assertEqual(container_image_path(image, prefix).read_bytes(), b"SIF")

    def test_get_envs_uses_environment_prefix_and_singularity_fallback(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "images"

            def find_runtime(executable):
                return "/usr/bin/singularity" if executable == "singularity" else None

            with (
                patch("caddsv.cli.shutil.which", side_effect=find_runtime),
                patch(
                    "caddsv.cli.subprocess.run",
                    side_effect=self.write_fake_image,
                ) as run_mock,
            ):
                result = self.runner.invoke(
                    app,
                    ["get", "envs", "--coordinate-based-only"],
                    env={"CADD_SV_SINGULARITY_PREFIX": str(prefix)},
                )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertTrue(prefix.is_dir())
            self.assertTrue(
                all(
                    call.args[0][0] == "/usr/bin/singularity"
                    for call in run_mock.call_args_list
                )
            )

    def test_get_envs_stops_and_cleans_up_after_pull_failure(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            prefix = Path(tmpdir) / "images"
            call_count = 0

            def fail_second_pull(command, check):
                nonlocal call_count
                call_count += 1
                if call_count == 2:
                    raise subprocess.CalledProcessError(1, command)
                self.write_fake_image(command, check)

            with (
                patch("caddsv.cli.shutil.which", return_value="/usr/bin/apptainer"),
                patch("caddsv.cli.subprocess.run", side_effect=fail_second_pull),
            ):
                result = self.runner.invoke(
                    app,
                    ["get", "envs", "--apptainer-prefix", str(prefix)],
                )

            self.assertEqual(result.exit_code, 2, result.output)
            self.assertEqual(call_count, 2)
            self.assertIn("Failed to download sv image", result.output)
            self.assertTrue(
                container_image_path(
                    DEFAULT_CONTAINER_IMAGES["preprocessing"], prefix
                ).is_file()
            )
            self.assertFalse(list(prefix.glob(".caddsv-image-*")))

    def test_get_envs_requires_a_container_runtime(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            result = self.runner.invoke(
                app,
                ["get", "envs", "--apptainer-prefix", tmpdir],
                env={"PATH": ""},
            )

        self.assertEqual(result.exit_code, 2, result.output)
        self.assertIn("Downloading environment images requires Apptainer", result.output)
        self.assertIn("Singularity.", result.output)


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
            source = workflow_output / "samplebed_score100.bed"
            source_contents = (
                "chr\tstart\tend\ttype\tCADD-SV_PHRED\tCADD-SV_score\n"
                "chr1\t1\t2\tDEL\t12.345678\t0.123456\n"
            )
            source.write_text(source_contents)

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
            self.assertNotEqual(
                run_mock.call_args.kwargs["env"].get("PYTHONNOUSERSITE"),
                "1",
            )
            self.assertFalse(cache_home.exists())
            self.assertEqual(source.read_text(), source_contents)
            self.assertEqual(
                (results / "scored" / "sample_score.tsv").read_text(),
                "#chr\tstart\tend\ttype\tCADD-SV_PHRED\tCADD-SV_score\n"
                "chr1\t1\t2\tDEL\t12.345678\t0.1235\n",
            )

    def test_sequence_only_mode_formats_final_scores_without_changing_header(self):
        runner = CliRunner()

        with tempfile.TemporaryDirectory() as tmpdir:
            root = Path(tmpdir)
            results = root / "results"
            sequence_input = root / "sample.tsv"
            sequence_input.write_text("A\tT\tDEL\n")

            workflow_output = results / "beds" / "sample" / "output"
            workflow_output.mkdir(parents=True)
            (workflow_output / "sample_seqonly_score.tsv").write_text(
                "id\ttype\tCADD-SV_seqonly_PHRED\tCADD-SV_seqonly_score\n"
                "sample\tDEL\t8.987654\t0.5\n"
            )

            with patch("caddsv.cli.subprocess.run"):
                result = runner.invoke(
                    app,
                    [
                        "run",
                        str(sequence_input),
                        "--seqonly",
                        "--output-dir",
                        str(results),
                        "--no-use-conda",
                    ],
                )

            self.assertEqual(result.exit_code, 0, result.output)
            self.assertEqual(
                (results / "scored" / "sample_seqonly_score.tsv").read_text(),
                "id\ttype\tCADD-SV_seqonly_PHRED\tCADD-SV_seqonly_score\n"
                "sample\tDEL\t8.987654\t0.5000\n",
            )

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
