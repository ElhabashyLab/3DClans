"""Unit tests for clans3d.utils.structure_utils."""
import os
import pytest
import pandas as pd
import requests
from unittest.mock import MagicMock, patch, call
from Bio.PDB.MMCIFParser import MMCIFParser

from clans3d.core.input_file_type import InputFileType
from clans3d.utils.structure_utils import (
    _log_interval,
    _parse_region_from_tsv_row,
    _fetch_alphafold_cif_url,
    download_alphafold_structure,
    extract_region_of_protein,
    fetch_structures,
)

MODULE = "clans3d.utils.structure_utils"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _tsv_row(start, end) -> pd.Series:
    """Build a pd.Series mimicking a TSV row with the given region values."""
    return pd.Series({"entry": "P00000", "region_start": start, "region_end": end})


def _tsv_row_no_region() -> pd.Series:
    return pd.Series({"entry": "P00000", "region_start": pd.NA, "region_end": pd.NA})


def _mock_response(json_data, status_ok=True):
    """Build a mock requests.Response."""
    resp = MagicMock()
    resp.json.return_value = json_data
    if not status_ok:
        resp.raise_for_status.side_effect = requests.HTTPError("400 Bad Request")
    return resp


# ---------------------------------------------------------------------------
# _log_interval
# ---------------------------------------------------------------------------

class TestLogInterval:
    def test_zero_total_returns_one(self):
        assert _log_interval(0) == 1

    def test_single_item_returns_one(self):
        assert _log_interval(1) == 1

    def test_small_total_clamps_to_one(self):
        # 4 // 5 == 0, so max(1, ...) → 1
        assert _log_interval(4) == 1

    def test_medium_total(self):
        # 50 // 5 == 10, min(10, 10) == 10
        assert _log_interval(50) == 10

    def test_large_total_capped_at_ten(self):
        # 1000 // 5 == 200, min(10, 200) == 10
        assert _log_interval(1000) == 10

    def test_25_items(self):
        # 25 // 5 == 5, min(10, 5) == 5
        assert _log_interval(25) == 5


# ---------------------------------------------------------------------------
# _parse_region_from_tsv_row
# ---------------------------------------------------------------------------

class TestParseRegionFromTsvRow:
    def test_both_missing_returns_none(self):
        assert _parse_region_from_tsv_row(_tsv_row_no_region()) is None

    def test_valid_region_returns_tuple(self):
        assert _parse_region_from_tsv_row(_tsv_row(1, 50)) == (1, 50)

    def test_valid_region_large_numbers(self):
        assert _parse_region_from_tsv_row(_tsv_row(100, 500)) == (100, 500)

    def test_only_start_missing_raises(self):
        row = pd.Series({"entry": "P00000", "region_start": pd.NA, "region_end": 50})
        with pytest.raises(ValueError, match="only one"):
            _parse_region_from_tsv_row(row)

    def test_only_end_missing_raises(self):
        row = pd.Series({"entry": "P00000", "region_start": 1, "region_end": pd.NA})
        with pytest.raises(ValueError, match="only one"):
            _parse_region_from_tsv_row(row)

    def test_start_zero_raises(self):
        with pytest.raises(ValueError, match="invalid region"):
            _parse_region_from_tsv_row(_tsv_row(0, 50))

    def test_start_negative_raises(self):
        with pytest.raises(ValueError, match="invalid region"):
            _parse_region_from_tsv_row(_tsv_row(-1, 50))

    def test_end_zero_raises(self):
        with pytest.raises(ValueError, match="invalid region"):
            _parse_region_from_tsv_row(_tsv_row(1, 0))

    def test_end_less_than_start_raises(self):
        with pytest.raises(ValueError, match="invalid region"):
            _parse_region_from_tsv_row(_tsv_row(50, 10))

    def test_start_equals_end_is_valid(self):
        # Single-residue region is biologically valid
        assert _parse_region_from_tsv_row(_tsv_row(5, 5)) == (5, 5)


# ---------------------------------------------------------------------------
# _fetch_alphafold_cif_url
# ---------------------------------------------------------------------------

class TestFetchAlphafoldCifUrl:
    CIF_URL = "https://alphafold.ebi.ac.uk/files/AF-P00000-F1-model_v4.cif"

    def test_returns_url_on_success(self):
        resp = _mock_response([{"cifUrl": self.CIF_URL}])
        with patch(f"{MODULE}.requests.get", return_value=resp):
            result = _fetch_alphafold_cif_url("P00000")
        assert result == self.CIF_URL

    def test_returns_none_on_http_error(self):
        resp = _mock_response([], status_ok=False)
        with patch(f"{MODULE}.requests.get", return_value=resp):
            result = _fetch_alphafold_cif_url("P00000")
        assert result is None

    def test_returns_none_on_connection_error(self):
        with patch(f"{MODULE}.requests.get", side_effect=requests.ConnectionError()):
            result = _fetch_alphafold_cif_url("P00000")
        assert result is None

    def test_returns_none_on_empty_predictions(self):
        resp = _mock_response([])
        with patch(f"{MODULE}.requests.get", return_value=resp):
            result = _fetch_alphafold_cif_url("P00000")
        assert result is None

    def test_returns_none_when_no_cif_url_in_response(self):
        resp = _mock_response([{"pdbUrl": "https://example.com/something.pdb"}])
        with patch(f"{MODULE}.requests.get", return_value=resp):
            result = _fetch_alphafold_cif_url("P00000")
        assert result is None

    def test_queries_correct_api_url(self):
        resp = _mock_response([{"cifUrl": self.CIF_URL}])
        with patch(f"{MODULE}.requests.get", return_value=resp) as mock_get:
            _fetch_alphafold_cif_url("P12345")
        mock_get.assert_called_once_with(
            "https://alphafold.ebi.ac.uk/api/prediction/P12345",
            timeout=30,
        )


# ---------------------------------------------------------------------------
# download_alphafold_structure
# ---------------------------------------------------------------------------

class TestDownloadAlphafoldStructure:
    CIF_URL = "https://alphafold.ebi.ac.uk/files/AF-P00000-F1-model_v4.cif"

    def test_returns_false_when_no_cif_url(self, tmp_path):
        with patch(f"{MODULE}._fetch_alphafold_cif_url", return_value=None):
            result = download_alphafold_structure("P00000", str(tmp_path), None)
        assert result is False

    def test_returns_false_when_download_fails(self, tmp_path):
        with patch(f"{MODULE}._fetch_alphafold_cif_url", return_value=self.CIF_URL), \
             patch(f"{MODULE}.download_file", return_value=False):
            result = download_alphafold_structure("P00000", str(tmp_path), None)
        assert result is False

    def test_returns_true_without_region(self, tmp_path):
        with patch(f"{MODULE}._fetch_alphafold_cif_url", return_value=self.CIF_URL), \
             patch(f"{MODULE}.download_file", return_value=True), \
             patch(f"{MODULE}.extract_region_of_protein") as mock_extract:
            result = download_alphafold_structure("P00000", str(tmp_path), None)
        assert result is True
        mock_extract.assert_not_called()

    def test_returns_true_with_region_and_calls_extract(self, tmp_path):
        region = (10, 50)
        expected_save_path = os.path.join(str(tmp_path), "P00000.cif")
        with patch(f"{MODULE}._fetch_alphafold_cif_url", return_value=self.CIF_URL), \
             patch(f"{MODULE}.download_file", return_value=True), \
             patch(f"{MODULE}.extract_region_of_protein") as mock_extract:
            result = download_alphafold_structure("P00000", str(tmp_path), region)
        assert result is True
        mock_extract.assert_called_once_with(expected_save_path, region, expected_save_path)

    def test_saves_to_correct_path(self, tmp_path):
        captured = {}
        def capture_download(url, path):
            captured["path"] = path
            return True
        with patch(f"{MODULE}._fetch_alphafold_cif_url", return_value=self.CIF_URL), \
             patch(f"{MODULE}.download_file", side_effect=capture_download), \
             patch(f"{MODULE}.extract_region_of_protein"):
            download_alphafold_structure("P12345", str(tmp_path), None)
        assert captured["path"] == os.path.join(str(tmp_path), "P12345.cif")


# ---------------------------------------------------------------------------
# fetch_structures (dispatch logic)
# ---------------------------------------------------------------------------

class TestFetchStructures:
    EXPECTED_RESULT = {"P00000": (1, 50)}

    def test_dispatches_fasta_to_process_fasta_file(self, small_fasta_path, tmp_path):
        with patch(f"{MODULE}.reset_dir_content"), \
             patch(f"{MODULE}.process_fasta_file", return_value=self.EXPECTED_RESULT) as mock_proc:
            result = fetch_structures(small_fasta_path, InputFileType.FASTA, str(tmp_path))
        mock_proc.assert_called_once_with(small_fasta_path, str(tmp_path))
        assert result == self.EXPECTED_RESULT

    def test_dispatches_a2m_to_process_fasta_file(self, small_fasta_path, tmp_path):
        with patch(f"{MODULE}.reset_dir_content"), \
             patch(f"{MODULE}.process_fasta_file", return_value=self.EXPECTED_RESULT) as mock_proc:
            result = fetch_structures(small_fasta_path, InputFileType.A2M, str(tmp_path))
        mock_proc.assert_called_once_with(small_fasta_path, str(tmp_path))
        assert result == self.EXPECTED_RESULT

    def test_dispatches_a3m_to_process_fasta_file(self, small_fasta_path, tmp_path):
        with patch(f"{MODULE}.reset_dir_content"), \
             patch(f"{MODULE}.process_fasta_file", return_value=self.EXPECTED_RESULT) as mock_proc:
            result = fetch_structures(small_fasta_path, InputFileType.A3M, str(tmp_path))
        mock_proc.assert_called_once_with(small_fasta_path, str(tmp_path))
        assert result == self.EXPECTED_RESULT

    def test_dispatches_tsv_to_process_tsv_file(self, small_tsv_path, tmp_path):
        with patch(f"{MODULE}.reset_dir_content"), \
             patch(f"{MODULE}.process_tsv_file", return_value=self.EXPECTED_RESULT) as mock_proc:
            result = fetch_structures(small_tsv_path, InputFileType.TSV, str(tmp_path))
        mock_proc.assert_called_once_with(small_tsv_path, str(tmp_path))
        assert result == self.EXPECTED_RESULT

    def test_raises_for_unsupported_type(self, tmp_path):
        with patch(f"{MODULE}.reset_dir_content"):
            with pytest.raises(ValueError, match="Unsupported"):
                fetch_structures("any.clans", InputFileType.CLANS, str(tmp_path))

    def test_resets_output_dir_before_processing(self, small_fasta_path, tmp_path):
        with patch(f"{MODULE}.reset_dir_content") as mock_reset, \
             patch(f"{MODULE}.process_fasta_file", return_value={}):
            fetch_structures(small_fasta_path, InputFileType.FASTA, str(tmp_path))
        mock_reset.assert_called_once_with(str(tmp_path))


# ---------------------------------------------------------------------------
# extract_region_of_protein
# ---------------------------------------------------------------------------

class TestExtractRegionOfProtein:
    def _count_residues(self, cif_path: str) -> int:
        parser = MMCIFParser(QUIET=True)
        structure = parser.get_structure("out", cif_path)
        return sum(1 for _ in structure.get_residues())

    def test_extracts_correct_number_of_residues(self, small_cif_path, tmp_path):
        out = str(tmp_path / "out.cif")
        # small.cif has residues 1-5; extract 2-4 → 3 residues
        extract_region_of_protein(small_cif_path, (2, 4), out)
        assert self._count_residues(out) == 3

    def test_full_range_keeps_all_residues(self, small_cif_path, tmp_path):
        out = str(tmp_path / "out.cif")
        extract_region_of_protein(small_cif_path, (1, 5), out)
        assert self._count_residues(out) == 5

    def test_single_residue_region(self, small_cif_path, tmp_path):
        out = str(tmp_path / "out.cif")
        extract_region_of_protein(small_cif_path, (3, 3), out)
        assert self._count_residues(out) == 1

    def test_default_output_path_has_roi_suffix(self, small_cif_path):
        result = extract_region_of_protein(small_cif_path, (1, 3))
        root, _ = os.path.splitext(small_cif_path)
        assert result == f"{root}_roi.cif"

    def test_custom_output_path_is_used(self, small_cif_path, tmp_path):
        out = str(tmp_path / "custom_out.cif")
        result = extract_region_of_protein(small_cif_path, (1, 3), out)
        assert result == out
        assert os.path.exists(out)
