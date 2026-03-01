"""Unit tests for clans3d.utils.api_utils._chunked."""
import pytest
from clans3d.utils.api_utils import _chunked


class TestChunked:
    def test_even_split(self):
        result = list(_chunked(["a", "b", "c", "d"], 2))
        assert result == [["a", "b"], ["c", "d"]]

    def test_uneven_split(self):
        result = list(_chunked(["a", "b", "c"], 2))
        assert result == [["a", "b"], ["c"]]

    def test_empty_list(self):
        result = list(_chunked([], 3))
        assert result == []

    def test_chunk_larger_than_list(self):
        result = list(_chunked(["x", "y"], 10))
        assert result == [["x", "y"]]

    def test_chunk_size_one(self):
        result = list(_chunked([1, 2, 3], 1))
        assert result == [[1], [2], [3]]

    def test_single_element(self):
        result = list(_chunked(["only"], 5))
        assert result == [["only"]]
