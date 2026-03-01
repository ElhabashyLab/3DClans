"""Unit tests for clans3d.utils.file_utils."""
import os
import pytest
from unittest.mock import patch, MagicMock
from clans3d.utils.file_utils import download_file, reset_dir_content, copy_dir_content


class TestResetDirContent:
    def test_creates_directory_if_missing(self, tmp_path):
        target = tmp_path / "new_dir"
        assert not target.exists()
        reset_dir_content(str(target))
        assert target.exists()
        assert target.is_dir()

    def test_removes_existing_files(self, tmp_path):
        target = tmp_path / "existing"
        target.mkdir()
        (target / "file.txt").write_text("data")
        reset_dir_content(str(target))
        assert target.exists()
        assert list(target.iterdir()) == []

    def test_removes_nested_subdirectory(self, tmp_path):
        target = tmp_path / "dir"
        target.mkdir()
        nested = target / "sub"
        nested.mkdir()
        (nested / "f.txt").write_text("x")
        reset_dir_content(str(target))
        assert target.exists()
        assert list(target.iterdir()) == []


class TestCopyDirContent:
    def test_copies_files(self, tmp_path):
        src = tmp_path / "src"
        src.mkdir()
        (src / "a.txt").write_text("hello")
        dst = tmp_path / "dst"
        copy_dir_content(str(src), str(dst))
        assert (dst / "a.txt").read_text() == "hello"

    def test_creates_dst_if_missing(self, tmp_path):
        src = tmp_path / "src"
        src.mkdir()
        (src / "a.txt").write_text("hi")
        dst = tmp_path / "nonexistent_dst"
        copy_dir_content(str(src), str(dst))
        assert dst.exists()

    def test_copies_nested_directories(self, tmp_path):
        src = tmp_path / "src"
        sub = src / "nested"
        sub.mkdir(parents=True)
        (sub / "b.txt").write_text("nested")
        dst = tmp_path / "dst"
        copy_dir_content(str(src), str(dst))
        assert (dst / "nested" / "b.txt").read_text() == "nested"


class TestDownloadFile:
    def test_returns_true_and_writes_file_on_success(self, tmp_path):
        out_path = str(tmp_path / "output" / "file.cif")
        mock_response = MagicMock()
        mock_response.raise_for_status = MagicMock()
        mock_response.iter_content = MagicMock(return_value=[b"chunk1", b"chunk2"])
        with patch("clans3d.utils.file_utils.requests.get", return_value=mock_response):
            result = download_file("http://example.com/file.cif", out_path)
        assert result is True
        assert os.path.exists(out_path)
        with open(out_path, "rb") as f:
            assert f.read() == b"chunk1chunk2"

    def test_returns_false_on_http_error(self, tmp_path):
        out_path = str(tmp_path / "file.cif")
        mock_response = MagicMock()
        mock_response.raise_for_status.side_effect = Exception("HTTP 404")
        with patch("clans3d.utils.file_utils.requests.get", return_value=mock_response):
            result = download_file("http://example.com/bad", out_path)
        assert result is False

    def test_returns_false_on_connection_error(self, tmp_path):
        out_path = str(tmp_path / "file.cif")
        with patch("clans3d.utils.file_utils.requests.get", side_effect=ConnectionError("refused")):
            result = download_file("http://example.com/file", out_path)
        assert result is False
