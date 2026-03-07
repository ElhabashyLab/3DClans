"""
Regression tests for CLANS file generation output stability.

These tests lock in the structural properties of generated CLANS files so
that refactors cannot silently break output shape. A fixed input FASTA and
fixed mock scores are used — no external binaries or network calls.
"""

import os

import pandas as pd
import pytest

from clans3d.core.clans_file_generator import ClansFileGenerator

FIXTURES_DIR = os.path.join(os.path.dirname(__file__), "..", "fixtures")
SMALL_FASTA = os.path.abspath(os.path.join(FIXTURES_DIR, "small.fasta"))

FIXED_SCORES = pd.DataFrame(
    {
        "Sequence_ID_1": ["P11111", "P11111", "P22222"],
        "Sequence_ID_2": ["P22222", "P33333", "P33333"],
        "score": [1.5e-30, 2.0e-28, 3.5e-25],
    }
)


@pytest.fixture
def generated_clans(tmp_path) -> str:
    """Generate a CLANS file from fixed inputs and return its path."""
    generator = ClansFileGenerator()
    out = str(tmp_path / "out.clans")
    generator.generate_clans_file(FIXED_SCORES, SMALL_FASTA, out)
    return out


class TestClansOutputStructure:
    def test_sequence_count_stable(self, generated_clans):
        content = open(generated_clans).read()
        assert "sequences=3" in content

    def test_coordinate_count_stable(self, generated_clans):
        content = open(generated_clans).read()
        pos_section = content.split("<pos>")[1].split("</pos>")[0]
        coord_lines = [l for l in pos_section.strip().splitlines() if l.strip()]
        assert len(coord_lines) == 3

    def test_score_entry_count_stable(self, generated_clans):
        content = open(generated_clans).read()
        hsp_section = content.split("<hsp>")[1].split("</hsp>")[0]
        score_lines = [l for l in hsp_section.strip().splitlines() if l.strip()]
        assert len(score_lines) == 3

    def test_all_required_sections_present(self, generated_clans):
        content = open(generated_clans).read()
        for section in (
            "<param>", "</param>",
            "<seq>", "</seq>",
            "<pos>", "</pos>",
            "<hsp>", "</hsp>",
        ):
            assert section in content

    def test_score_format_stable(self, generated_clans):
        """Each score line must follow the 'idx1 idx2:score' format."""
        import re
        content = open(generated_clans).read()
        hsp_section = content.split("<hsp>")[1].split("</hsp>")[0]
        score_lines = [l.strip() for l in hsp_section.strip().splitlines() if l.strip()]
        pattern = re.compile(r"^\d+ \d+:[0-9eE+\-.]+$")
        for line in score_lines:
            assert pattern.match(line), f"Unexpected score line format: {line!r}"

    def test_fasta_sequences_preserved(self, generated_clans):
        """All three sequence records from the input FASTA must appear in <seq>."""
        content = open(generated_clans).read()
        seq_section = content.split("<seq>")[1].split("</seq>")[0]
        for uid in ("P11111", "P22222", "P33333"):
            assert uid in seq_section
