"""Unit tests for clans3d.similarity.foldseek.Foldseek."""
import subprocess
import pytest
import pandas as pd
from unittest.mock import patch
from clans3d.similarity.foldseek import Foldseek
from clans3d.similarity.tm_mode import TmMode


@pytest.fixture
def foldseek_evalue():
    return Foldseek(score="evalue", tm_mode=TmMode.MIN)


@pytest.fixture
def foldseek_tm():
    return Foldseek(score="TM", tm_mode=TmMode.MIN)


@pytest.fixture
def foldseek_tm_max():
    return Foldseek(score="TM", tm_mode=TmMode.MAX)


@pytest.fixture
def foldseek_tm_avg():
    return Foldseek(score="TM", tm_mode=TmMode.AVG)


# ---------------------------------------------------------------------------
# _get_output_columns_from_self_score
# ---------------------------------------------------------------------------

class TestGetOutputColumns:
    def test_evalue_columns(self, foldseek_evalue):
        assert foldseek_evalue.arg_output_columns == "query,target,evalue"

    def test_tm_columns(self, foldseek_tm):
        assert foldseek_tm.arg_output_columns == "query,target,qtmscore,ttmscore"

    def test_invalid_score_raises(self):
        with pytest.raises(ValueError):
            Foldseek(score="rmsd", tm_mode=TmMode.MIN)


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
    def test_score_is_one_minus_min_tm_by_default(self, foldseek_tm):
        df = pd.DataFrame({
            "Sequence_ID_1": ["A"],
            "Sequence_ID_2": ["B"],
            "TM1": [0.7],
            "TM2": [0.9],
        })
        result = foldseek_tm._clean_scores(df)
        assert result["score"].iloc[0] == pytest.approx(1.0 - 0.7)

    def test_score_is_one_minus_max_tm_when_mode_max(self, foldseek_tm_max):
        df = pd.DataFrame({
            "Sequence_ID_1": ["A"],
            "Sequence_ID_2": ["B"],
            "TM1": [0.7],
            "TM2": [0.9],
        })
        result = foldseek_tm_max._clean_scores(df)
        assert result["score"].iloc[0] == pytest.approx(1.0 - 0.9)

    def test_score_is_one_minus_avg_tm_when_mode_avg(self, foldseek_tm_avg):
        df = pd.DataFrame({
            "Sequence_ID_1": ["A"],
            "Sequence_ID_2": ["B"],
            "TM1": [0.7],
            "TM2": [0.9],
        })
        result = foldseek_tm_avg._clean_scores(df)
        assert result["score"].iloc[0] == pytest.approx(1.0 - 0.8)

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


# ---------------------------------------------------------------------------
# _create_database
# ---------------------------------------------------------------------------

class TestCreateDatabase:
    def test_raises_runtime_error_on_subprocess_failure(self, foldseek_evalue, tmp_path):
        foldseek_evalue.working_dir = str(tmp_path)
        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(1, "foldseek", stderr="db error")
            with pytest.raises(RuntimeError, match="createdb failed"):
                foldseek_evalue._create_database(str(tmp_path / "structures"), "testDB")

    def test_returns_correct_db_path_on_success(self, foldseek_evalue, tmp_path):
        foldseek_evalue.working_dir = str(tmp_path)
        with patch("subprocess.run"):
            result = foldseek_evalue._create_database(str(tmp_path / "structures"), "testDB")
        assert result == str(tmp_path / "testDB")


# ---------------------------------------------------------------------------
# _parse_output
# ---------------------------------------------------------------------------

class TestParseOutputEvalue:
    def test_raises_on_convertalis_subprocess_failure(self, foldseek_evalue, tmp_path):
        foldseek_evalue.working_dir = str(tmp_path)
        foldseek_evalue.db = str(tmp_path / "db")
        foldseek_evalue.alignmentDb = str(tmp_path / "aln")
        with patch("subprocess.run") as mock_run:
            mock_run.side_effect = subprocess.CalledProcessError(1, "foldseek", stderr="conv error")
            with pytest.raises(RuntimeError, match="convertalis failed"):
                foldseek_evalue._parse_output()

    def test_parses_evalue_output_correctly(self, foldseek_evalue, tmp_path):
        foldseek_evalue.working_dir = str(tmp_path)
        foldseek_evalue.db = str(tmp_path / "db")
        foldseek_evalue.alignmentDb = str(tmp_path / "aln")
        # pandas uses row 0 as the column header, so write a dummy header row
        # followed by one real data row (matches real foldseek convertalis output)
        (tmp_path / "alignmentFile").write_text("A\tB\t1e-5\nC\tD\t1e-10\n")
        with patch("subprocess.run"):
            result = foldseek_evalue._parse_output()
        assert list(result.columns) == ["Sequence_ID_1", "Sequence_ID_2", "score"]
        assert result["score"].iloc[0] == pytest.approx(1e-10)


class TestParseOutputTM:
    def test_parses_tm_output_correctly(self, foldseek_tm, tmp_path):
        foldseek_tm.working_dir = str(tmp_path)
        foldseek_tm.db = str(tmp_path / "db")
        foldseek_tm.alignmentDb = str(tmp_path / "aln")
        # pandas uses row 0 as the column header; row 1 is the actual data
        (tmp_path / "alignmentFile").write_text("A\tB\t0.5\t0.6\nC\tD\t0.7\t0.9\n")
        with patch("subprocess.run"):
            result = foldseek_tm._parse_output()
        # default mode is min: min(0.7, 0.9) = 0.7 -> score = 1 - 0.7 = 0.3
        assert result["score"].iloc[0] == pytest.approx(0.3)
        assert "TM1" not in result.columns
        assert "TM2" not in result.columns
