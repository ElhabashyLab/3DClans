"""Unit tests for clans3d.similarity.foldseek.Foldseek."""
import pytest
import pandas as pd
from clans3d.similarity.foldseek import Foldseek


@pytest.fixture
def foldseek_evalue():
    return Foldseek(score="evalue")


@pytest.fixture
def foldseek_tm():
    return Foldseek(score="TM")


# ---------------------------------------------------------------------------
# _get_output_columns_from_self_score
# ---------------------------------------------------------------------------

class TestGetOutputColumns:
    def test_evalue_columns(self, foldseek_evalue):
        assert foldseek_evalue.output_columns == "query,target,evalue"

    def test_tm_columns(self, foldseek_tm):
        assert foldseek_tm.output_columns == "query,target,qtmscore,ttmscore"

    def test_invalid_score_raises(self):
        with pytest.raises(ValueError):
            Foldseek(score="rmsd")


# ---------------------------------------------------------------------------
# _clean_scores  (evalue mode)
# ---------------------------------------------------------------------------

class TestCleanScoresEvalue:
    def test_self_hits_removed(self, foldseek_evalue):
        df = pd.DataFrame({
            "Sequence_ID_1": ["A", "A", "B"],
            "Sequence_ID_2": ["A", "B", "B"],  # A-A and B-B are self-hits
            "evalue": [0.0, 1e-10, 0.0],
        })
        result = foldseek_evalue._clean_scores(df)
        assert all(result["Sequence_ID_1"] != result["Sequence_ID_2"])

    def test_symmetric_duplicates_removed(self, foldseek_evalue):
        df = pd.DataFrame({
            "Sequence_ID_1": ["A", "B"],
            "Sequence_ID_2": ["B", "A"],
            "evalue": [1e-10, 1e-10],
        })
        result = foldseek_evalue._clean_scores(df)
        assert len(result) == 1

    def test_score_column_equals_evalue(self, foldseek_evalue):
        df = pd.DataFrame({
            "Sequence_ID_1": ["A"],
            "Sequence_ID_2": ["B"],
            "evalue": [1e-20],
        })
        result = foldseek_evalue._clean_scores(df)
        assert "score" in result.columns
        assert result["score"].iloc[0] == pytest.approx(1e-20)

    def test_no_evalue_column_in_output(self, foldseek_evalue):
        df = pd.DataFrame({
            "Sequence_ID_1": ["A"],
            "Sequence_ID_2": ["B"],
            "evalue": [1e-5],
        })
        result = foldseek_evalue._clean_scores(df)
        assert "evalue" not in result.columns


# ---------------------------------------------------------------------------
# _clean_scores  (TM mode)
# ---------------------------------------------------------------------------

class TestCleanScoresTM:
    def test_score_is_one_minus_max_tm(self, foldseek_tm):
        df = pd.DataFrame({
            "Sequence_ID_1": ["A"],
            "Sequence_ID_2": ["B"],
            "TM1": [0.7],
            "TM2": [0.9],
        })
        result = foldseek_tm._clean_scores(df)
        assert result["score"].iloc[0] == pytest.approx(1.0 - 0.9)

    def test_tm_columns_removed_from_output(self, foldseek_tm):
        df = pd.DataFrame({
            "Sequence_ID_1": ["A"],
            "Sequence_ID_2": ["B"],
            "TM1": [0.5],
            "TM2": [0.6],
        })
        result = foldseek_tm._clean_scores(df)
        assert "TM1" not in result.columns
        assert "TM2" not in result.columns

    def test_self_hits_removed(self, foldseek_tm):
        df = pd.DataFrame({
            "Sequence_ID_1": ["A", "A"],
            "Sequence_ID_2": ["A", "B"],
            "TM1": [1.0, 0.8],
            "TM2": [1.0, 0.7],
        })
        result = foldseek_tm._clean_scores(df)
        assert len(result) == 1


# ---------------------------------------------------------------------------
# _detect_phase
# ---------------------------------------------------------------------------

class TestDetectPhase:
    def test_detects_prefiltering(self, foldseek_evalue):
        reader = {"chunks": ["some line", "Running prefiltering step..."]}
        assert foldseek_evalue._detect_phase(reader) == "Prefiltering ..."

    def test_detects_structurealign(self, foldseek_evalue):
        reader = {"chunks": ["StructureAlign running structurealign"]}
        assert foldseek_evalue._detect_phase(reader) == "Structure alignment ..."

    def test_returns_empty_when_no_match(self, foldseek_evalue):
        reader = {"chunks": ["no relevant output here"]}
        assert foldseek_evalue._detect_phase(reader) == ""

    def test_returns_latest_phase_from_reversed_chunks(self, foldseek_evalue):
        # convertalis appears later → should be returned
        reader = {"chunks": ["counting k-mers step", "convertalis running"]}
        phase = foldseek_evalue._detect_phase(reader)
        assert phase == "Converting results ..."
