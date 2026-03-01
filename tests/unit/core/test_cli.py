"""Unit tests for clans3d.core.cli.parse_args."""
import os
import pytest
from clans3d.core.cli import parse_args


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
