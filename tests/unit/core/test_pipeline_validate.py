"""Unit tests for clans3d.core.pipeline.ClansPipeline._validate, PipelineConfig defaults, and fetch_structures."""
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

    def test_default_output_filename_is_none(self, fasta_path):
        config = PipelineConfig(
            input_file=fasta_path,
            input_type=InputFileType.FASTA,
            tool=ToolType.FOLDSEEK,
        )
        assert config.output_filename is None

    def test_custom_output_filename(self, fasta_path):
        config = PipelineConfig(
            input_file=fasta_path,
            input_type=InputFileType.FASTA,
            tool=ToolType.FOLDSEEK,
            output_filename="custom.clans",
        )
        assert config.output_filename == "custom.clans"


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


class TestFetchStructures:
    @pytest.fixture
    def pipeline(self, fasta_path):
        config = PipelineConfig(
            input_file=fasta_path,
            input_type=InputFileType.FASTA,
            tool=ToolType.FOLDSEEK,
        )
        return ClansPipeline(config)

    def test_raises_when_no_structures_downloaded(self, pipeline):
        """Empty dict from structure_utils must abort with RuntimeError."""
        with patch("clans3d.core.pipeline.fetch_structures", return_value={}):
            with pytest.raises(RuntimeError, match="No structures"):
                pipeline.fetch_structures()

    def test_returns_uid_dict_on_partial_success(self, pipeline):
        """A non-empty dict is returned unchanged — partial failures are fine."""
        expected = {"P11111": None, "P22222": (1, 100)}
        with patch("clans3d.core.pipeline.fetch_structures", return_value=expected), \
             patch("os.makedirs"):
            result = pipeline.fetch_structures()
        assert result == expected


class TestGenerateCleanedFastaTSV:
    def test_forwards_download_workers_to_tsv_fasta_generation(self, tsv_path, tmp_path):
        config = PipelineConfig(
            input_file=tsv_path,
            input_type=InputFileType.TSV,
            tool=ToolType.FOLDSEEK,
            cleaned_input_storage=str(tmp_path / "cleaned"),
            download_workers=7,
        )
        pipeline = ClansPipeline(config)

        uids_with_regions = {"P11111": None}
        expected_path = str(tmp_path / "cleaned" / "input_cleaned.fasta")

        with patch("os.makedirs"), \
             patch("clans3d.core.pipeline.generate_fasta_from_uids_with_regions", return_value=expected_path) as mock_generate:
            result = pipeline.generate_cleaned_fasta(uids_with_regions)

        assert result == expected_path
        mock_generate.assert_called_once_with(
            uids_with_regions,
            expected_path,
            max_workers=7,
        )
