"""Unit tests for clans3d.core.clans_file.ClansFile."""
import pytest
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from clans3d.core.clans_file import ClansFile


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_records():
    return [
        SeqRecord(Seq("ACDE"), id="sp|P11111|PROT1_HUMAN/1-50", description="sp|P11111|PROT1_HUMAN/1-50"),
        SeqRecord(Seq("FGHI"), id="sp|P22222|PROT2_MOUSE/10-60", description="sp|P22222|PROT2_MOUSE/10-60"),
        SeqRecord(Seq("KLMN"), id="sp|P33333|PROT3_RAT", description="sp|P33333|PROT3_RAT"),
    ]


def _make_scores():
    return pd.DataFrame({
        "Sequence_ID_1": ["P11111", "P11111", "P22222"],
        "Sequence_ID_2": ["P22222", "P33333", "P33333"],
        "score": [1e-30, 2e-28, 3e-25],
    })


def _make_coords():
    return [
        ("P11111", 0.1, 0.2, 0.3),
        ("P22222", 0.4, 0.5, 0.6),
        ("P33333", 0.7, 0.8, 0.9),
    ]


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

class TestClansFileConstruction:
    def test_constructs_from_fasta_records(self):
        cf = ClansFile(3, _make_coords(), _make_scores(), fasta_records=_make_records())
        assert cf.number_of_sequences == 3

    def test_raises_when_neither_records_nor_path_given(self):
        with pytest.raises(ValueError):
            ClansFile(3, _make_coords(), _make_scores(),
                      fasta_records=None, path_to_fasta=None)

    def test_uid_to_index_built_correctly(self):
        cf = ClansFile(3, _make_coords(), _make_scores(), fasta_records=_make_records())
        assert cf._uid_to_index["P11111"] == 0
        assert cf._uid_to_index["P22222"] == 1
        assert cf._uid_to_index["P33333"] == 2


# ---------------------------------------------------------------------------
# _validate_data_consistency
# ---------------------------------------------------------------------------

class TestValidation:
    def test_raises_when_coord_count_wrong(self):
        coords_short = [("P11111", 0.1, 0.2, 0.3)]  # only 1 of 3
        with pytest.raises(ValueError, match="coordinates"):
            ClansFile(3, coords_short, _make_scores(), fasta_records=_make_records())

    def test_raises_when_coord_uid_not_in_records(self):
        bad_coords = [
            ("UNKNOWN", 0.1, 0.2, 0.3),
            ("P22222",  0.4, 0.5, 0.6),
            ("P33333",  0.7, 0.8, 0.9),
        ]
        with pytest.raises(ValueError, match="UNKNOWN"):
            ClansFile(3, bad_coords, _make_scores(), fasta_records=_make_records())

    def test_raises_when_score_uid_not_in_records(self):
        bad_scores = pd.DataFrame({
            "Sequence_ID_1": ["P11111"],
            "Sequence_ID_2": ["GHOST"],
            "score": [1e-10],
        })
        with pytest.raises(ValueError, match="GHOST"):
            ClansFile(3, _make_coords(), bad_scores, fasta_records=_make_records())

    def test_empty_scores_df_accepted(self):
        empty = pd.DataFrame(columns=["Sequence_ID_1", "Sequence_ID_2", "score"])
        cf = ClansFile(3, _make_coords(), empty, fasta_records=_make_records())
        assert cf.scores.empty


# ---------------------------------------------------------------------------
# __str__
# ---------------------------------------------------------------------------

class TestClansFileStr:
    def setup_method(self):
        self.cf = ClansFile(3, _make_coords(), _make_scores(), fasta_records=_make_records())
        self.content = str(self.cf)

    def test_contains_sequence_count(self):
        assert "sequences=3" in self.content

    def test_all_sections_present(self):
        for section in ["<param>", "</param>", "<seq>", "</seq>", "<pos>", "</pos>", "<hsp>", "</hsp>"]:
            assert section in self.content

    def test_coordinates_use_indices(self):
        # P11111 is index 0
        assert "0 0.1 0.2 0.3" in self.content

    def test_scores_use_indices(self):
        # P11111=0, P22222=1
        assert "0 1:" in self.content

    def test_fasta_sequences_present(self):
        assert "ACDE" in self.content
        assert "FGHI" in self.content

    def test_parameters_written_with_dash_prefix(self):
        cf = ClansFile(3, _make_coords(), _make_scores(), fasta_records=_make_records(),
                       parameters={"pval": "1e-5", "maxmove": "0.1"})
        content = str(cf)
        assert "-pval 1e-5" in content
        assert "-maxmove 0.1" in content

    def test_empty_scores_produces_empty_hsp_block(self):
        empty = pd.DataFrame(columns=["Sequence_ID_1", "Sequence_ID_2", "score"])
        cf = ClansFile(3, _make_coords(), empty, fasta_records=_make_records())
        content = str(cf)
        assert "<hsp>" in content
        assert "</hsp>" in content
        between = content.split("<hsp>")[1].split("</hsp>")[0].strip()
        assert between == ""
