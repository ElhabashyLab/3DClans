"""Unit tests for clans3d.similarity.usalign.USalign."""
import pytest
import pandas as pd
from clans3d.similarity.usalign import USalign


USALIGN_OUTPUT = (
    "#PDBchain1\tPDBchain2\tTM1\tTM2\tRMSD\tID1\tID2\tIDali\tL1\tL2\tLali\n"
    "P11111.cif\tP22222.cif\t0.75\t0.80\t1.2\t0.5\t0.5\t0.5\t100\t110\t90\n"
    "P11111.cif\tP33333.cif\t0.60\t0.65\t2.1\t0.4\t0.4\t0.4\t100\t120\t85\n"
    "P22222.cif\tP33333.cif\t0.50\t0.55\t2.5\t0.3\t0.3\t0.3\t110\t120\t80\n"
)


@pytest.fixture
def usalign():
    return USalign()


class TestCountDataLines:
    def test_ignores_header_lines(self, usalign):
        reader = {"chunks": ["#Header line\n", "P1\tP2\t0.8\n"]}
        count = usalign._count_data_lines(reader)
        assert count == 1

    def test_ignores_warning_lines(self, usalign):
        reader = {"chunks": ["Warning: something\n", "P1\tP2\t0.8\n"]}
        count = usalign._count_data_lines(reader)
        assert count == 1

    def test_counts_multiple_data_lines(self, usalign):
        reader = {"chunks": ["P1\tP2\t0.8\n", "P1\tP3\t0.7\n", "P2\tP3\t0.6\n"]}
        count = usalign._count_data_lines(reader)
        assert count == 3

    def test_ignores_blank_lines(self, usalign):
        reader = {"chunks": ["\n", "  \n", "P1\tP2\t0.8\n"]}
        count = usalign._count_data_lines(reader)
        assert count == 1


class TestParseOutput:
    def test_returns_dataframe_with_correct_columns(self, usalign):
        usalign.output = USALIGN_OUTPUT
        result = usalign._parse_output()
        assert set(["Sequence_ID_1", "Sequence_ID_2", "score", "TM"]).issubset(result.columns)

    def test_extensions_stripped_from_ids(self, usalign):
        usalign.output = USALIGN_OUTPUT
        result = usalign._parse_output()
        assert all("." not in uid for uid in result["Sequence_ID_1"])
        assert all("." not in uid for uid in result["Sequence_ID_2"])

    def test_score_is_one_minus_max_tm(self, usalign):
        usalign.output = USALIGN_OUTPUT
        result = usalign._parse_output()
        row = result[result["Sequence_ID_1"] == "P11111"].iloc[0]
        # TM1=0.75, TM2=0.80 → max=0.80 → score=0.20
        assert row["score"] == pytest.approx(1.0 - 0.80)

    def test_correct_number_of_rows(self, usalign):
        usalign.output = USALIGN_OUTPUT
        result = usalign._parse_output()
        assert len(result) == 3

    def test_empty_output_raises(self, usalign):
        usalign.output = ""
        with pytest.raises(RuntimeError):
            usalign._parse_output()

    def test_tm_column_not_in_final_output(self, usalign):
        """Raw TM1/TM2 columns should be dropped; the intermediate TM column remains for reference."""
        usalign.output = USALIGN_OUTPUT
        result = usalign._parse_output()
        assert "TM1" not in result.columns
        assert "TM2" not in result.columns

    def test_raises_on_malformed_output(self, usalign):
        """Non-TSV / completely garbled output should raise RuntimeError."""
        usalign.output = "this is not\tvalid usalign\noutput at all"
        with pytest.raises(RuntimeError):
            usalign._parse_output()
