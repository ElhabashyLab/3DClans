"""Unit tests for clans3d.core.pipeline.ClansPipeline._validate and PipelineConfig defaults."""
import os
import pytest
from unittest.mock import patch
from clans3d.core.pipeline import ClansPipeline, PipelineConfig
from clans3d.core.input_file_type import InputFileType
from clans3d.similarity.tool_type import ToolType


@pytest.fixture
def fasta_path(tmp_path):
    p = tmp_path / "input.fasta"
    p.write_text(">P11111\nACDE\n")
    return str(p)


@pytest.fixture
def tsv_path(tmp_path):
    p = tmp_path / "input.tsv"
    p.write_text("entry\tregion_start\tregion_end\nP11111\t1\t50\n")
    return str(p)


class TestPipelineConfigDefaults:
    def test_default_structures_dir(self, fasta_path):
        config = PipelineConfig(
            input_file=fasta_path,
            input_type=InputFileType.FASTA,
            tool=ToolType.FOLDSEEK,
        )
        assert "structures" in config.structures_dir

    def test_default_output_dir(self, fasta_path):
        config = PipelineConfig(
            input_file=fasta_path,
            input_type=InputFileType.FASTA,
            tool=ToolType.FOLDSEEK,
        )
        assert "clans_files" in config.output_dir


class TestClansPipelineValidate:
    def test_raises_file_not_found(self, tmp_path):
        config = PipelineConfig(
            input_file=str(tmp_path / "nonexistent.fasta"),
            input_type=InputFileType.FASTA,
            tool=ToolType.FOLDSEEK,
        )
        with pytest.raises(FileNotFoundError):
            ClansPipeline(config)

    def test_raises_on_extension_mismatch(self, tsv_path):
        config = PipelineConfig(
            input_file=tsv_path,
            input_type=InputFileType.FASTA,  # mismatch: file is .tsv
            tool=ToolType.FOLDSEEK,
        )
        with pytest.raises(ValueError):
            ClansPipeline(config)

    def test_raises_when_foldseek_score_set_for_usalign(self, fasta_path):
        config = PipelineConfig(
            input_file=fasta_path,
            input_type=InputFileType.FASTA,
            tool=ToolType.USALIGN,
            foldseek_score="evalue",
        )
        with pytest.raises(ValueError, match="foldseek_score"):
            ClansPipeline(config)

    def test_valid_config_does_not_raise(self, fasta_path):
        config = PipelineConfig(
            input_file=fasta_path,
            input_type=InputFileType.FASTA,
            tool=ToolType.FOLDSEEK,
        )
        pipeline = ClansPipeline(config)  # should not raise
        assert pipeline is not None
