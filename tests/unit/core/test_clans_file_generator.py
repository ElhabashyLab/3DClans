"""Unit tests for clans3d.core.clans_file_generator.ClansFileGenerator."""
import io
import pytest
import pandas as pd
from Bio import SeqIO
from clans3d.core.clans_file_generator import ClansFileGenerator


MINIMAL_CLANS_LINES = [
    "sequences=3",
    "<param>",
    "pval=1e-5",
    "</param>",
    "<seq>",
    ">sp|P11111|PROT1_HUMAN/1-50",
    "ACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLM",
    ">sp|P22222|PROT2_MOUSE/10-60",
    "WYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIK",
    ">sp|P33333|PROT3_RAT",
    "MNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYAC",
    "</seq>",
    "<pos>",
    "0 0.100 0.200 0.300",
    "1 0.400 0.500 0.600",
    "2 0.700 0.800 0.900",
    "</pos>",
    "<hsp>",
    "0 1:1.5e-30",
    "0 2:2.0e-28",
    "1 2:3.5e-25",
    "</hsp>",
]


@pytest.fixture
def generator():
    return ClansFileGenerator()


# ---------------------------------------------------------------------------
# _extract_block
# ---------------------------------------------------------------------------

class TestExtractBlock:
    def test_extracts_correct_lines(self, generator):
        block = generator._extract_block("param", MINIMAL_CLANS_LINES)
        assert block == ["pval=1e-5"]

    def test_returns_none_for_missing_tag(self, generator):
        result = generator._extract_block("missing", MINIMAL_CLANS_LINES)
        assert result is None


# ---------------------------------------------------------------------------
# _parse_number_of_sequences
# ---------------------------------------------------------------------------

class TestParseNumberOfSequences:
    def test_parses_correctly(self, generator):
        assert generator._parse_number_of_sequences(MINIMAL_CLANS_LINES) == 3

    def test_raises_when_missing(self, generator):
        with pytest.raises(ValueError, match="sequences="):
            generator._parse_number_of_sequences(["<param>", "</param>"])


# ---------------------------------------------------------------------------
# _parse_param_block
# ---------------------------------------------------------------------------

class TestParseParamBlock:
    def test_returns_dict(self, generator):
        params = generator._parse_param_block(MINIMAL_CLANS_LINES)
        assert params == {"pval": "1e-5"}

    def test_returns_none_when_missing(self, generator):
        lines_without_param = [l for l in MINIMAL_CLANS_LINES if "param" not in l and "pval" not in l]
        result = generator._parse_param_block(lines_without_param)
        assert result is None


# ---------------------------------------------------------------------------
# _parse_fasta_block
# ---------------------------------------------------------------------------

class TestParseFastaBlock:
    def test_returns_three_records(self, generator):
        records = generator._parse_fasta_block(MINIMAL_CLANS_LINES)
        assert len(records) == 3

    def test_raises_when_seq_block_missing(self, generator):
        no_seq = [l for l in MINIMAL_CLANS_LINES if "<seq>" not in l and "</seq>" not in l
                  and not l.startswith(">") and "ACDE" not in l
                  and "WYAC" not in l and "MNPQ" not in l]
        with pytest.raises(ValueError):
            generator._parse_fasta_block(no_seq)


# ---------------------------------------------------------------------------
# _parse_pos_block
# ---------------------------------------------------------------------------

class TestParsePosBlock:
    def test_returns_tuples(self, generator):
        positions = generator._parse_pos_block(MINIMAL_CLANS_LINES)
        assert positions[0] == (0, 0.1, 0.2, 0.3)
        assert len(positions) == 3

    def test_raises_when_missing(self, generator):
        no_pos = [l for l in MINIMAL_CLANS_LINES if "pos" not in l and "0.1" not in l
                  and "0.4" not in l and "0.7" not in l]
        with pytest.raises(ValueError):
            generator._parse_pos_block(no_pos)


# ---------------------------------------------------------------------------
# _parse_scores_block
# ---------------------------------------------------------------------------

class TestParseScoresBlock:
    def test_returns_dataframe_with_correct_columns(self, generator):
        df = generator._parse_scores_block(MINIMAL_CLANS_LINES)
        assert list(df.columns) == ["Sequence_ID_1", "Sequence_ID_2", "score"]

    def test_correct_number_of_rows(self, generator):
        df = generator._parse_scores_block(MINIMAL_CLANS_LINES)
        assert len(df) == 3

    def test_score_values_are_floats(self, generator):
        df = generator._parse_scores_block(MINIMAL_CLANS_LINES)
        assert df["score"].dtype == float

    def test_raises_when_hsp_block_missing(self, generator):
        no_hsp = [l for l in MINIMAL_CLANS_LINES if "hsp" not in l
                  and "1:1.5" not in l and "2:2.0" not in l and "2:3.5" not in l]
        with pytest.raises(ValueError):
            generator._parse_scores_block(no_hsp)


# ---------------------------------------------------------------------------
# _normalize_scores_format
# ---------------------------------------------------------------------------

class TestNormalizeScoresFormat:
    def test_removes_symmetric_duplicates(self, generator):
        df = pd.DataFrame({
            "Sequence_ID_1": ["A", "B"],
            "Sequence_ID_2": ["B", "A"],
            "score": [0.5, 0.5],
        })
        result = generator._normalize_scores_format(df, ["A", "B"])
        assert len(result) == 1

    def test_preserves_asymmetric_pairs(self, generator):
        df = pd.DataFrame({
            "Sequence_ID_1": ["A", "A"],
            "Sequence_ID_2": ["B", "C"],
            "score": [0.5, 0.3],
        })
        result = generator._normalize_scores_format(df, ["A", "B", "C"])
        assert len(result) == 2


# ---------------------------------------------------------------------------
# _generate_random_coordinates
# ---------------------------------------------------------------------------

class TestGenerateRandomCoordinates:
    def test_returns_one_tuple_per_uid(self, generator):
        uids = ["A", "B", "C"]
        coords = generator._generate_random_coordinates(uids)
        assert len(coords) == 3

    def test_all_xyz_in_zero_one(self, generator):
        coords = generator._generate_random_coordinates(["X", "Y", "Z"])
        for uid, x, y, z in coords:
            assert 0.0 <= x <= 1.0
            assert 0.0 <= y <= 1.0
            assert 0.0 <= z <= 1.0

    def test_uids_preserved(self, generator):
        uids = ["P1", "P2"]
        coords = generator._generate_random_coordinates(uids)
        returned_uids = [c[0] for c in coords]
        assert returned_uids == uids


# ---------------------------------------------------------------------------
# parse_clans_file (end-to-end parsing from fixture file)
# ---------------------------------------------------------------------------

class TestParseClansFile:
    def test_parses_fixture_correctly(self, generator, small_clans_path):
        cf = generator.parse_clans_file(small_clans_path)
        assert cf.number_of_sequences == 3
        assert len(cf.fasta_records) == 3
        assert len(cf.coordinates) == 3
        assert len(cf.scores) == 3

    def test_uid_mapping_correct(self, generator, small_clans_path):
        cf = generator.parse_clans_file(small_clans_path)
        coord_uids = {c[0] for c in cf.coordinates}
        assert "P11111" in coord_uids
        assert "P22222" in coord_uids
        assert "P33333" in coord_uids
