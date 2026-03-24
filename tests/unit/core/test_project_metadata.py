"""Unit tests for project metadata in pyproject.toml."""
import tomllib
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[3]
PYPROJECT_PATH = REPO_ROOT / "pyproject.toml"


def test_project_name_is_3dclans():
    data = tomllib.loads(PYPROJECT_PATH.read_text())
    assert data["project"]["name"] == "3dclans"


def test_console_script_name_is_3dclans():
    data = tomllib.loads(PYPROJECT_PATH.read_text())
    assert data["project"]["scripts"]["3dclans"] == "clans3d.main:main"
