"""Unit tests for clans3d.core.cli.parse_args and resolve_output_path."""
import os
import pytest
from clans3d.core.cli import parse_args, resolve_output_path


@pytest.fixture
def fasta_path(tmp_path):
    p = tmp_path / "proteins.fasta"
    p.write_text(">P1\nACDE\n")
    return str(p)


class TestParseArgs:
    def test_minimal_valid_args(self, fasta_path):
        args = parse_args(["-l", fasta_path, "-i", "fasta", "-t", "foldseek"])
        assert args.load == fasta_path
        assert args.input_type == "fasta"
        assert args.tool == "foldseek"

    def test_usalign_tool(self, fasta_path):
        args = parse_args(["-l", fasta_path, "-i", "fasta", "-t", "USalign"])
        assert args.tool == "USalign"

    def test_verbose_flag(self, fasta_path):
        args = parse_args(["-l", fasta_path, "-i", "fasta", "-t", "foldseek", "-v"])
        assert args.verbose is True
        assert args.quiet is False

    def test_quiet_flag(self, fasta_path):
        args = parse_args(["-l", fasta_path, "-i", "fasta", "-t", "foldseek", "-q"])
        assert args.quiet is True
        assert args.verbose is False

    def test_score_argument(self, fasta_path):
        args = parse_args(["-l", fasta_path, "-i", "fasta", "-t", "foldseek", "-s", "TM"])
        assert args.score == "TM"

    def test_default_score_is_none(self, fasta_path):
        args = parse_args(["-l", fasta_path, "-i", "fasta", "-t", "foldseek"])
        assert args.score is None

    def test_default_tm_mode_is_min(self, fasta_path):
        args = parse_args(["-l", fasta_path, "-i", "fasta", "-t", "foldseek"])
        assert args.tm_mode == "min"

    def test_tm_mode_argument(self, fasta_path):
        args = parse_args([
            "-l", fasta_path, "-i", "fasta", "-t", "foldseek", "--tm_mode", "mean"
        ])
        assert args.tm_mode == "mean"

    def test_invalid_tool_raises(self, fasta_path):
        with pytest.raises(SystemExit):
            parse_args(["-l", fasta_path, "-i", "fasta", "-t", "blast"])

    def test_invalid_input_type_raises(self, fasta_path):
        with pytest.raises(SystemExit):
            parse_args(["-l", fasta_path, "-i", "xml", "-t", "foldseek"])

    def test_config_file_merged(self, fasta_path, tmp_path):
        conf = tmp_path / "test.conf"
        conf.write_text(f"-load {fasta_path}\n-input_type fasta\n-tool foldseek\n")
        args = parse_args(["-l", fasta_path, "-i", "fasta", "-t", "foldseek", "-c", str(conf)])
        assert args.tool == "foldseek"

    def test_cli_overrides_config(self, fasta_path, tmp_path):
        conf = tmp_path / "test.conf"
        # config says evalue but CLI says TM
        conf.write_text("-score evalue\n")
        args = parse_args([
            "-l", fasta_path, "-i", "fasta", "-t", "foldseek",
            "-s", "TM", "-c", str(conf),
        ])
        assert args.score == "TM"

    def test_cli_tm_mode_overrides_config(self, fasta_path, tmp_path):
        conf = tmp_path / "test.conf"
        conf.write_text("-tm_mode min\n")
        args = parse_args([
            "-l", fasta_path, "-i", "fasta", "-t", "foldseek",
            "--tm_mode", "max", "-c", str(conf),
        ])
        assert args.tm_mode == "max"

    def test_default_out_is_none(self, fasta_path):
        args = parse_args(["-l", fasta_path, "-i", "fasta", "-t", "foldseek"])
        assert args.out is None

    def test_out_flag_with_file_path(self, fasta_path):
        args = parse_args([
            "-l", fasta_path, "-i", "fasta", "-t", "foldseek",
            "-o", "/tmp/results/my_output.clans",
        ])
        assert args.out == "/tmp/results/my_output.clans"

    def test_out_flag_with_directory(self, fasta_path):
        args = parse_args([
            "-l", fasta_path, "-i", "fasta", "-t", "foldseek",
            "-o", "/tmp/results/",
        ])
        assert args.out == "/tmp/results/"

    def test_out_flag_from_config_file(self, fasta_path, tmp_path):
        conf = tmp_path / "test.conf"
        conf.write_text("-out /tmp/from_config.clans\n")
        args = parse_args([
            "-l", fasta_path, "-i", "fasta", "-t", "foldseek",
            "-c", str(conf),
        ])
        assert args.out == "/tmp/from_config.clans"

    def test_cli_out_overrides_config_out(self, fasta_path, tmp_path):
        conf = tmp_path / "test.conf"
        conf.write_text("-out /tmp/from_config.clans\n")
        args = parse_args([
            "-l", fasta_path, "-i", "fasta", "-t", "foldseek",
            "-o", "/tmp/from_cli.clans", "-c", str(conf),
        ])
        assert args.out == "/tmp/from_cli.clans"


class TestResolveOutputPath:
    def test_none_returns_default(self):
        output_dir, output_filename = resolve_output_path(None)
        assert output_dir == os.path.join("output", "clans_files")
        assert output_filename is None

    def test_clans_file_path(self):
        output_dir, output_filename = resolve_output_path("/tmp/results/my_output.clans")
        assert output_dir == "/tmp/results"
        assert output_filename == "my_output.clans"

    def test_clans_file_without_dir(self):
        output_dir, output_filename = resolve_output_path("my_output.clans")
        assert output_dir == "."
        assert output_filename == "my_output.clans"

    def test_directory_path(self):
        output_dir, output_filename = resolve_output_path("/tmp/results/")
        assert output_dir == "/tmp/results/"
        assert output_filename is None

    def test_directory_without_trailing_slash(self):
        output_dir, output_filename = resolve_output_path("/tmp/results")
        assert output_dir == "/tmp/results"
        assert output_filename is None
