"""
Integration tests for ClansPipeline.run() with mocked structure fetching and
similarity computation.

External dependencies (AlphaFold downloads, Foldseek/USalign subprocesses) are
patched out so that only the pipeline orchestration logic and the CLANS file
generation are exercised with real I/O.
"""

import os

import pandas as pd
import pytest
from unittest.mock import patch
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from clans3d.core.pipeline import ClansPipeline, PipelineConfig
from clans3d.core.input_file_type import InputFileType
from clans3d.similarity.tool_type import ToolType


def _make_seq_record(uid: str, region=None) -> SeqRecord:
    """Create a minimal SeqRecord for mocking download_fasta_record."""
    seq = "ACDEFGHIKLMNPQRSTVWY" * 3
    record_id = uid if region is None else f"{uid}/{region[0]}-{region[1]}"
    return SeqRecord(Seq(seq), id=record_id, description=record_id)

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "..", "fixtures")
SMALL_FASTA = os.path.abspath(os.path.join(FIXTURES_DIR, "small.fasta"))
SMALL_TSV = os.path.abspath(os.path.join(FIXTURES_DIR, "small.tsv"))

# UIDs that match the three records in fixtures/small.fasta and fixtures/small.tsv
MOCK_UIDS_WITH_REGIONS = {
    "P11111": None,
    "P22222": None,
    "P33333": None,
}

MOCK_SCORES = pd.DataFrame(
    {
        "Sequence_ID_1": ["P11111", "P11111", "P22222"],
        "Sequence_ID_2": ["P22222", "P33333", "P33333"],
        "score": [1.5e-30, 2.0e-28, 3.5e-25],
    }
)


@pytest.fixture
def pipeline(tmp_path):
    config = PipelineConfig(
        input_file=SMALL_FASTA,
        input_type=InputFileType.FASTA,
        tool=ToolType.FOLDSEEK,
        structures_dir=str(tmp_path / "structures"),
        output_dir=str(tmp_path / "output"),
        cleaned_input_storage=str(tmp_path / "cleaned"),
    )
    return ClansPipeline(config)


@pytest.fixture
def tsv_pipeline(tmp_path):
    config = PipelineConfig(
        input_file=SMALL_TSV,
        input_type=InputFileType.TSV,
        tool=ToolType.FOLDSEEK,
        structures_dir=str(tmp_path / "structures"),
        output_dir=str(tmp_path / "output"),
        cleaned_input_storage=str(tmp_path / "cleaned"),
    )
    return ClansPipeline(config)


class TestPipelineMocked:
    def _run_mocked(self, pipeline):
        """Run the pipeline with fetch_structures and StructSimComputer mocked out."""
        with (
            patch(
                "clans3d.core.pipeline.fetch_structures",
                return_value=MOCK_UIDS_WITH_REGIONS,
            ),
            patch("clans3d.core.pipeline.StructSimComputer") as MockComputer,
        ):
            MockComputer.return_value.run.return_value = MOCK_SCORES
            return pipeline.run()

    def test_run_writes_clans_file(self, pipeline):
        clans_path, _ = self._run_mocked(pipeline)
        assert os.path.isfile(clans_path)

    def test_run_output_has_correct_sequence_count(self, pipeline):
        clans_path, _ = self._run_mocked(pipeline)
        content = open(clans_path).read()
        assert "sequences=3" in content

    def test_run_output_contains_all_hsp_entries(self, pipeline):
        clans_path, _ = self._run_mocked(pipeline)
        content = open(clans_path).read()
        hsp_section = content.split("<hsp>")[1].split("</hsp>")[0]
        score_lines = [l for l in hsp_section.strip().splitlines() if l.strip()]
        assert len(score_lines) == 3

    def test_run_returns_existing_fasta_path(self, pipeline):
        _, fasta_path = self._run_mocked(pipeline)
        assert os.path.isfile(fasta_path)
        assert fasta_path.endswith(".fasta")

    def test_run_output_contains_all_required_sections(self, pipeline):
        clans_path, _ = self._run_mocked(pipeline)
        content = open(clans_path).read()
        for section in ("<param>", "</param>", "<seq>", "</seq>", "<pos>", "</pos>", "<hsp>", "</hsp>"):
            assert section in content

    def test_run_raises_when_no_structures_available(self, pipeline):
        with (
            patch(
                "clans3d.core.pipeline.fetch_structures",
                return_value={},
            ),
        ):
            with pytest.raises(RuntimeError, match="No structures"):
                pipeline.run()

    def test_run_with_custom_output_filename(self, tmp_path):
        config = PipelineConfig(
            input_file=SMALL_FASTA,
            input_type=InputFileType.FASTA,
            tool=ToolType.FOLDSEEK,
            structures_dir=str(tmp_path / "structures"),
            output_dir=str(tmp_path / "custom_output"),
            cleaned_input_storage=str(tmp_path / "cleaned"),
            output_filename="my_results.clans",
        )
        custom_pipeline = ClansPipeline(config)
        with (
            patch(
                "clans3d.core.pipeline.fetch_structures",
                return_value=MOCK_UIDS_WITH_REGIONS,
            ),
            patch("clans3d.core.pipeline.StructSimComputer") as MockComputer,
        ):
            MockComputer.return_value.run.return_value = MOCK_SCORES
            clans_path, _ = custom_pipeline.run()
        assert os.path.isfile(clans_path)
        assert clans_path == os.path.join(str(tmp_path / "custom_output"), "my_results.clans")


class TestPipelineMockedTSV:
    """Integration tests for the TSV input path.

    For TSV inputs the pipeline does NOT copy sequences from the original file;
    instead it fetches them from UniProt via ``generate_fasta_from_uids_with_regions``.
    That network call is patched here so the test remains offline.
    """

    def _run_mocked(self, tsv_pipeline):
        with (
            patch(
                "clans3d.core.pipeline.fetch_structures",
                return_value=MOCK_UIDS_WITH_REGIONS,
            ),
            patch("clans3d.core.pipeline.StructSimComputer") as MockComputer,
            patch(
                "clans3d.utils.fasta_utils.uniprot_accessions_to_uniparc_accessions",
                return_value={"P11111": "UPI000", "P22222": "UPI001", "P33333": "UPI002"},
            ),
            patch(
                "clans3d.utils.fasta_utils.download_fasta_record",
                side_effect=lambda uid, upi, region: _make_seq_record(uid, region),
            ),
        ):
            MockComputer.return_value.run.return_value = MOCK_SCORES
            return tsv_pipeline.run()

    def test_run_writes_clans_file(self, tsv_pipeline):
        clans_path, _ = self._run_mocked(tsv_pipeline)
        assert os.path.isfile(clans_path)

    def test_run_output_has_correct_sequence_count(self, tsv_pipeline):
        clans_path, _ = self._run_mocked(tsv_pipeline)
        content = open(clans_path).read()
        assert "sequences=3" in content

    def test_run_output_contains_all_hsp_entries(self, tsv_pipeline):
        clans_path, _ = self._run_mocked(tsv_pipeline)
        content = open(clans_path).read()
        hsp_section = content.split("<hsp>")[1].split("</hsp>")[0]
        score_lines = [l for l in hsp_section.strip().splitlines() if l.strip()]
        assert len(score_lines) == 3

    def test_run_returns_existing_fasta_path(self, tsv_pipeline):
        _, fasta_path = self._run_mocked(tsv_pipeline)
        assert os.path.isfile(fasta_path)
        assert fasta_path.endswith(".fasta")


class TestPipelineMockedLocalStructuresDb:
    def test_run_uses_local_structures_db_branch(self, pipeline, tmp_path):
        structures_db = tmp_path / "structures_db"
        structures_db.mkdir()
        config = PipelineConfig(
            input_file=SMALL_FASTA,
            input_type=InputFileType.FASTA,
            tool=ToolType.FOLDSEEK,
            structures_dir=str(tmp_path / "structures"),
            output_dir=str(tmp_path / "output"),
            cleaned_input_storage=str(tmp_path / "cleaned"),
            structures_db=str(structures_db),
        )
        local_pipeline = ClansPipeline(config)
        with (
            patch(
                "clans3d.core.pipeline.prepare_structures_from_local_db",
                return_value=MOCK_UIDS_WITH_REGIONS,
            ) as mock_prepare,
            patch("clans3d.core.pipeline.fetch_structures") as mock_download,
            patch("clans3d.core.pipeline.StructSimComputer") as MockComputer,
        ):
            MockComputer.return_value.run.return_value = MOCK_SCORES
            clans_path, fasta_path = local_pipeline.run()

        assert os.path.isfile(clans_path)
        assert os.path.isfile(fasta_path)
        mock_prepare.assert_called_once()
        mock_download.assert_not_called()
