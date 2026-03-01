"""Unit tests for clans3d.similarity.struct_sim_computer.StructSimComputer."""
import pytest
import pandas as pd
from unittest.mock import MagicMock, patch
from clans3d.similarity.struct_sim_computer import StructSimComputer
from clans3d.similarity.tool_type import ToolType
from clans3d.similarity.foldseek import Foldseek
from clans3d.similarity.usalign import USalign


class TestCreateTool:
    def test_returns_foldseek_instance(self):
        computer = StructSimComputer(foldseek_score="evalue")
        tool = computer._create_tool(ToolType.FOLDSEEK)
        assert isinstance(tool, Foldseek)

    def test_returns_usalign_instance(self):
        computer = StructSimComputer()
        tool = computer._create_tool(ToolType.USALIGN)
        assert isinstance(tool, USalign)

    def test_raises_for_unknown_tool_type(self):
        computer = StructSimComputer()
        mock_type = MagicMock()
        mock_type.__eq__ = lambda self, other: False
        with pytest.raises(ValueError):
            computer._create_tool(mock_type)

    def test_foldseek_score_passed_to_tool(self):
        computer = StructSimComputer(foldseek_score="TM")
        tool = computer._create_tool(ToolType.FOLDSEEK)
        assert tool.score == "TM"


class TestRun:
    def test_logs_warning_on_score_count_mismatch(self, tmp_path):
        """When tool returns fewer scores than expected, a warning should be logged."""
        # n=4: expected = n*(n+1)//2 - n = 10-4 = 6
        struct_dir = tmp_path / "structures"
        struct_dir.mkdir()
        for i in range(4):
            (struct_dir / f"P{i}.cif").write_text("CIF data")

        # Only return 3 scores instead of 6
        short_df = pd.DataFrame({
            "Sequence_ID_1": ["A", "B", "C"],
            "Sequence_ID_2": ["B", "C", "D"],
            "score": [0.1, 0.2, 0.3],
        })
        computer = StructSimComputer(foldseek_score="evalue")
        with patch("clans3d.similarity.struct_sim_computer.logger") as mock_logger, \
             patch.object(computer, "_create_tool") as mock_factory:
            mock_tool = MagicMock()
            mock_tool.name = "MockTool"
            mock_tool.start_run.return_value = short_df
            mock_factory.return_value = mock_tool
            computer.run(ToolType.FOLDSEEK, str(struct_dir))
        mock_logger.warning.assert_called_once()

    def test_returns_dataframe_from_tool(self, tmp_path):
        struct_dir = tmp_path / "structures"
        struct_dir.mkdir()
        (struct_dir / "P1.cif").write_text("")
        (struct_dir / "P2.cif").write_text("")

        expected_df = pd.DataFrame({
            "Sequence_ID_1": ["P1"],
            "Sequence_ID_2": ["P2"],
            "score": [1e-10],
        })
        computer = StructSimComputer()
        with patch.object(computer, "_create_tool") as mock_factory:
            mock_tool = MagicMock()
            mock_tool.name = "MockTool"
            mock_tool.start_run.return_value = expected_df
            mock_factory.return_value = mock_tool
            result = computer.run(ToolType.FOLDSEEK, str(struct_dir))
        assert len(result) == 1
