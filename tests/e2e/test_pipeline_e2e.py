"""
End-to-end tests for the CLANS pipeline using real Foldseek and USalign binaries.

These tests download real AlphaFold structures, run the selected tool, and verify
that a valid CLANS file is produced. They are skipped automatically when the
required binary is not found in PATH.

Run with:
    pytest tests/e2e/ -m e2e
"""

import os
import re
import shutil
import subprocess

import pytest

from clans3d.core.pipeline import ClansPipeline, PipelineConfig
from clans3d.core.input_file_type import InputFileType
from clans3d.similarity.tool_type import ToolType

# 5-entry FASTA with 4 real UniProt accessions + 1 unavailable entry.
# The pipeline tolerates missing structures and produces output for those that succeed.
EXAMPLES_DIR = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "examples")
)
SMALL_FASTA = os.path.join(EXAMPLES_DIR, "small_fasta_files", "5.fasta")
SMALL_TSV = os.path.join(EXAMPLES_DIR, "small_tsv_files", "5.tsv")

# Expected number of successfully downloadable sequences (4 real MYOD1 orthologs)
MIN_EXPECTED_SEQUENCES = 3


def _parse_sequence_count(clans_path: str) -> int:
    with open(clans_path) as f:
        for line in f:
            m = re.match(r"sequences=(\d+)", line.strip())
            if m:
                return int(m.group(1))
    raise ValueError(f"No 'sequences=' line found in {clans_path}")


def _parse_score_count(clans_path: str) -> int:
    content = open(clans_path).read()
    if "<hsp>" not in content:
        return 0
    hsp_section = content.split("<hsp>")[1].split("</hsp>")[0]
    return len([l for l in hsp_section.strip().splitlines() if l.strip()])


@pytest.fixture
def foldseek_pipeline(tmp_path):
    if not shutil.which("foldseek"):
        pytest.skip("foldseek not in PATH")
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
def foldseek_tm_pipeline(tmp_path):
    if not shutil.which("foldseek"):
        pytest.skip("foldseek not in PATH")
    config = PipelineConfig(
        input_file=SMALL_FASTA,
        input_type=InputFileType.FASTA,
        tool=ToolType.FOLDSEEK,
        foldseek_score="TM",
        structures_dir=str(tmp_path / "structures"),
        output_dir=str(tmp_path / "output"),
        cleaned_input_storage=str(tmp_path / "cleaned"),
    )
    return ClansPipeline(config)


@pytest.fixture
def usalign_pipeline(tmp_path):
    if not shutil.which("USalign"):
        pytest.skip("USalign not in PATH")
    config = PipelineConfig(
        input_file=SMALL_FASTA,
        input_type=InputFileType.FASTA,
        tool=ToolType.USALIGN,
        structures_dir=str(tmp_path / "structures"),
        output_dir=str(tmp_path / "output"),
        cleaned_input_storage=str(tmp_path / "cleaned"),
    )
    return ClansPipeline(config)


@pytest.mark.e2e
class TestPipelineFoldseekEvalue:
    def test_clans_file_is_created(self, foldseek_pipeline):
        clans_path, _ = foldseek_pipeline.run()
        assert os.path.isfile(clans_path)

    def test_sequence_count_meets_minimum(self, foldseek_pipeline):
        clans_path, _ = foldseek_pipeline.run()
        assert _parse_sequence_count(clans_path) >= MIN_EXPECTED_SEQUENCES

    def test_score_entries_present(self, foldseek_pipeline):
        clans_path, _ = foldseek_pipeline.run()
        assert _parse_score_count(clans_path) >= 1

    def test_all_required_sections_present(self, foldseek_pipeline):
        clans_path, _ = foldseek_pipeline.run()
        content = open(clans_path).read()
        for section in ("<seq>", "</seq>", "<pos>", "</pos>", "<hsp>", "</hsp>"):
            assert section in content


@pytest.mark.e2e
class TestPipelineFoldseekTM:
    def test_clans_file_is_created(self, foldseek_tm_pipeline):
        clans_path, _ = foldseek_tm_pipeline.run()
        assert os.path.isfile(clans_path)

    def test_scores_are_in_unit_interval(self, foldseek_tm_pipeline):
        """TM-based scores are transformed to 1 - max(TM1, TM2) and must be in [0, 1]."""
        clans_path, _ = foldseek_tm_pipeline.run()
        content = open(clans_path).read()
        hsp_section = content.split("<hsp>")[1].split("</hsp>")[0]
        for line in hsp_section.strip().splitlines():
            line = line.strip()
            if not line:
                continue
            score = float(line.split(":")[1])
            assert 0.0 <= score <= 1.0, f"Score out of [0,1]: {score}"


@pytest.mark.e2e
class TestPipelineUSalign:
    def test_clans_file_is_created(self, usalign_pipeline):
        clans_path, _ = usalign_pipeline.run()
        assert os.path.isfile(clans_path)

    def test_sequence_count_meets_minimum(self, usalign_pipeline):
        clans_path, _ = usalign_pipeline.run()
        assert _parse_sequence_count(clans_path) >= MIN_EXPECTED_SEQUENCES

    def test_score_entries_present(self, usalign_pipeline):
        clans_path, _ = usalign_pipeline.run()
        assert _parse_score_count(clans_path) >= 1


@pytest.mark.e2e
class TestCliEntrypoint:
    def test_cli_exits_with_zero(self, tmp_path):
        if not shutil.which("foldseek"):
            pytest.skip("foldseek not in PATH")
        result = subprocess.run(
            [
                "3dclans",
                "-l", SMALL_FASTA,
                "-i", "fasta",
                "-t", "foldseek",
            ],
            capture_output=True,
            text=True,
            cwd=str(tmp_path),  # default output and work dirs go into tmp
        )
        assert result.returncode == 0, (
            f"3dclans exited with {result.returncode}\n"
            f"stdout: {result.stdout}\n"
            f"stderr: {result.stderr}"
        )


@pytest.fixture
def foldseek_tsv_pipeline(tmp_path):
    if not shutil.which("foldseek"):
        pytest.skip("foldseek not in PATH")
    config = PipelineConfig(
        input_file=SMALL_TSV,
        input_type=InputFileType.TSV,
        tool=ToolType.FOLDSEEK,
        structures_dir=str(tmp_path / "structures"),
        output_dir=str(tmp_path / "output"),
        cleaned_input_storage=str(tmp_path / "cleaned"),
    )
    return ClansPipeline(config)


@pytest.mark.e2e
class TestPipelineTSVInput:
    """E2E tests for the TSV input path.

    TSV input differs from FASTA: the cleaned FASTA is generated from scratch
    by downloading sequences from UniProt rather than copying from the original file.
    """

    def test_clans_file_is_created(self, foldseek_tsv_pipeline):
        clans_path, _ = foldseek_tsv_pipeline.run()
        assert os.path.isfile(clans_path)

    def test_sequence_count_meets_minimum(self, foldseek_tsv_pipeline):
        clans_path, _ = foldseek_tsv_pipeline.run()
        assert _parse_sequence_count(clans_path) >= MIN_EXPECTED_SEQUENCES

    def test_score_entries_present(self, foldseek_tsv_pipeline):
        clans_path, _ = foldseek_tsv_pipeline.run()
        assert _parse_score_count(clans_path) >= 1

    def test_cleaned_fasta_is_written(self, foldseek_tsv_pipeline):
        _, fasta_path = foldseek_tsv_pipeline.run()
        assert os.path.isfile(fasta_path)
