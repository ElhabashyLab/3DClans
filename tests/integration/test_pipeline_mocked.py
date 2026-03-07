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

from clans3d.core.pipeline import ClansPipeline, PipelineConfig
from clans3d.core.input_file_type import InputFileType
from clans3d.similarity.tool_type import ToolType

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "..", "fixtures")
SMALL_FASTA = os.path.abspath(os.path.join(FIXTURES_DIR, "small.fasta"))

# UIDs that match the three records in fixtures/small.fasta
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
