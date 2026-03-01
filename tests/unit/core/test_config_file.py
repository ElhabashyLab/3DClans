"""Unit tests for clans3d.core.config_file.ConfigFile."""
import os
import pytest
from clans3d.core.config_file import ConfigFile


class TestWriteAndReadConfig:
    def test_round_trip(self, tmp_path):
        path = str(tmp_path / "test.conf")
        cfg = ConfigFile(path)
        data = {"load": "proteins.fasta", "tool": "foldseek", "score": "evalue"}
        cfg.write_config(data)
        result = cfg.read_config()
        assert result == data

    def test_written_file_uses_dash_prefix(self, tmp_path):
        path = str(tmp_path / "test.conf")
        cfg = ConfigFile(path)
        cfg.write_config({"key": "value"})
        content = open(path).read()
        assert "-key value" in content

    def test_documentation_written_as_comment(self, tmp_path):
        path = str(tmp_path / "doc.conf")
        cfg = ConfigFile(path, documentation="test doc")
        cfg.write_config({"k": "v"})
        content = open(path).read()
        assert content.startswith("# test doc")

    def test_ignores_blank_lines_and_comments(self, tmp_path):
        path = str(tmp_path / "commented.conf")
        with open(path, "w") as f:
            f.write("# ignore this\n\n-tool foldseek\n")
        cfg = ConfigFile(path)
        result = cfg.read_config()
        assert result == {"tool": "foldseek"}

    def test_raises_on_missing_dash_prefix(self, tmp_path):
        path = str(tmp_path / "bad.conf")
        with open(path, "w") as f:
            f.write("tool foldseek\n")
        cfg = ConfigFile(path)
        with pytest.raises(ValueError):
            cfg.read_config()

    def test_raises_on_empty_key(self, tmp_path):
        path = str(tmp_path / "emptykey.conf")
        with open(path, "w") as f:
            f.write("- value\n")
        cfg = ConfigFile(path)
        with pytest.raises(ValueError):
            cfg.read_config()

    def test_empty_value_allowed(self, tmp_path):
        path = str(tmp_path / "novalue.conf")
        with open(path, "w") as f:
            f.write("-flag\n")
        cfg = ConfigFile(path)
        result = cfg.read_config()
        assert result == {"flag": ""}


class TestConfigToArgv:
    def test_returns_double_dash_list(self, tmp_path):
        path = str(tmp_path / "argv.conf")
        cfg = ConfigFile(path)
        cfg.write_config({"tool": "foldseek", "load": "file.fasta"})
        argv = cfg.config_to_argv()
        assert "--tool" in argv
        assert "foldseek" in argv
        assert "--load" in argv
        assert "file.fasta" in argv

    def test_pairs_are_consecutive(self, tmp_path):
        path = str(tmp_path / "argv.conf")
        cfg = ConfigFile(path)
        cfg.write_config({"key": "val"})
        argv = cfg.config_to_argv()
        idx = argv.index("--key")
        assert argv[idx + 1] == "val"
