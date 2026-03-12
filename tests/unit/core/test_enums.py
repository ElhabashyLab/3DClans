"""Unit tests for clans3d.core.input_file_type and clans3d.similarity.tool_type enums."""
import pytest
from clans3d.core.input_file_type import InputFileType
from clans3d.similarity.tool_type import ToolType


class TestInputFileType:
    @pytest.mark.parametrize("value, expected", [
        ("fasta", InputFileType.FASTA),
        ("a2m",   InputFileType.A2M),
        ("a3m",   InputFileType.A3M),
        ("clans", InputFileType.CLANS),
    ])
    def test_from_string(self, value, expected):
        assert InputFileType(value) == expected

    def test_invalid_value_raises(self):
        with pytest.raises(ValueError):
            InputFileType("xml")

    def test_value_attribute(self):
        assert InputFileType.FASTA.value == "fasta"


class TestToolType:
    @pytest.mark.parametrize("value, expected", [
        ("foldseek", ToolType.FOLDSEEK),
        ("USalign",  ToolType.USALIGN),
    ])
    def test_from_string(self, value, expected):
        assert ToolType(value) == expected

    def test_invalid_value_raises(self):
        with pytest.raises(ValueError):
            ToolType("blast")

    def test_value_attribute(self):
        assert ToolType.FOLDSEEK.value == "foldseek"
        assert ToolType.USALIGN.value == "USalign"
