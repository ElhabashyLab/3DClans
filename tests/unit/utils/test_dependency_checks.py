"""Unit tests for clans3d.utils.dependency_checks."""
import pytest
from unittest.mock import patch
from clans3d.utils.dependency_checks import check_external_tool, verify_tool_dependencies
from clans3d.similarity.tool_type import ToolType


class TestCheckExternalTool:
    def test_returns_path_when_found(self):
        with patch("shutil.which", return_value="/usr/bin/foldseek"):
            result = check_external_tool("foldseek")
        assert result == "/usr/bin/foldseek"

    def test_returns_none_when_not_found(self):
        with patch("shutil.which", return_value=None):
            result = check_external_tool("nonexistent_tool")
        assert result is None


class TestVerifyToolDependencies:
    def test_raises_runtime_error_when_missing(self):
        with patch("clans3d.utils.dependency_checks.check_external_tool", return_value=None):
            with pytest.raises(RuntimeError, match="foldseek"):
                verify_tool_dependencies(ToolType.FOLDSEEK)

    def test_no_error_when_foldseek_present(self):
        with patch("clans3d.utils.dependency_checks.check_external_tool", return_value="/usr/bin/foldseek"):
            verify_tool_dependencies(ToolType.FOLDSEEK)  # should not raise

    def test_raises_for_usalign_missing(self):
        with patch("clans3d.utils.dependency_checks.check_external_tool", return_value=None):
            with pytest.raises(RuntimeError, match="USalign"):
                verify_tool_dependencies(ToolType.USALIGN)
