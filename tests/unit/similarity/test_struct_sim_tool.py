"""Unit tests for clans3d.similarity.struct_sim_tool.StructSimTool._execute_run."""
import pytest
import pandas as pd
from clans3d.similarity.struct_sim_tool import StructSimTool


class _ConcreteSimTool(StructSimTool):
    """Minimal concrete subclass for testing base-class behaviour without external tools."""

    def __init__(self):
        super().__init__("test-tool", "A test tool", "work/test")
        self._parse_output_result = pd.DataFrame()

    def _parse_output(self) -> pd.DataFrame:
        return self._parse_output_result

    def _log_progress(self, stdout_reader: dict, stderr_reader: dict) -> None:
        pass


@pytest.fixture
def tool():
    return _ConcreteSimTool()


# ---------------------------------------------------------------------------
# _execute_run
# ---------------------------------------------------------------------------

class TestExecuteRun:
    def test_raises_runtime_error_on_nonzero_exit(self, tool):
        tool.command = ["sh", "-c", "exit 1"]
        with pytest.raises(RuntimeError, match="subprocess failed"):
            tool._execute_run(0)

    def test_output_captured_on_success(self, tool):
        tool.command = ["echo", "hello"]
        tool._execute_run(0)
        assert "hello" in tool.output

    def test_returns_dataframe_from_parse_output(self, tool):
        expected = pd.DataFrame({
            "Sequence_ID_1": ["A"],
            "Sequence_ID_2": ["B"],
            "score": [0.5],
        })
        tool._parse_output_result = expected
        tool.command = ["echo", ""]
        result = tool._execute_run(0)
        pd.testing.assert_frame_equal(result, expected)

    def test_expected_number_of_scores_stored(self, tool):
        tool.command = ["echo", ""]
        tool._execute_run(42)
        assert tool.expected_number_of_scores == 42

    def test_error_message_contains_exit_code(self, tool):
        tool.command = ["sh", "-c", "exit 2"]
        with pytest.raises(RuntimeError, match="exit 2"):
            tool._execute_run(0)
