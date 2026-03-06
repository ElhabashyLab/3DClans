"""Unit tests for clans3d.utils.api_utils."""
import pytest
import requests
from unittest.mock import patch, MagicMock, call

from clans3d.utils.api_utils import (
    _chunked,
    _submit_idmapping_job,
    _wait_for_results,
    _fetch_paginated_results,
)

MODULE = "clans3d.utils.api_utils"
RUN_URL = "https://rest.uniprot.org/idmapping/run"
RESULTS_URL = "https://rest.uniprot.org/idmapping/results/{}"


def _mock_response(json_data=None, status_code=200, links=None):
    resp = MagicMock()
    resp.status_code = status_code
    resp.json.return_value = json_data or {}
    resp.links = links or {}
    if status_code >= 400:
        resp.raise_for_status.side_effect = requests.HTTPError(f"{status_code}")
    else:
        resp.raise_for_status.return_value = None
    return resp


# ---------------------------------------------------------------------------
# _chunked
# ---------------------------------------------------------------------------

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


# ---------------------------------------------------------------------------
# _submit_idmapping_job
# ---------------------------------------------------------------------------

class TestSubmitIdmappingJob:
    def test_returns_job_id(self):
        resp = _mock_response({"jobId": "abc123"})
        with patch(f"{MODULE}.requests.post", return_value=resp):
            job_id = _submit_idmapping_job(["P11111", "P22222"], RUN_URL)
        assert job_id == "abc123"

    def test_posts_to_correct_url(self):
        resp = _mock_response({"jobId": "x"})
        with patch(f"{MODULE}.requests.post", return_value=resp) as mock_post:
            _submit_idmapping_job(["P11111"], RUN_URL)
        mock_post.assert_called_once()
        assert mock_post.call_args[0][0] == RUN_URL

    def test_sends_correct_from_to_fields(self):
        resp = _mock_response({"jobId": "x"})
        with patch(f"{MODULE}.requests.post", return_value=resp) as mock_post:
            _submit_idmapping_job(["P11111"], RUN_URL)
        data = mock_post.call_args[1]["data"]
        assert data["from"] == "UniProtKB_AC-ID"
        assert data["to"] == "UniParc"

    def test_joins_accessions_as_comma_separated(self):
        resp = _mock_response({"jobId": "x"})
        with patch(f"{MODULE}.requests.post", return_value=resp) as mock_post:
            _submit_idmapping_job(["P11111", "P22222", "P33333"], RUN_URL)
        data = mock_post.call_args[1]["data"]
        assert data["ids"] == "P11111,P22222,P33333"

    def test_raises_on_http_error(self):
        resp = _mock_response(status_code=400)
        with patch(f"{MODULE}.requests.post", return_value=resp):
            with pytest.raises(requests.HTTPError):
                _submit_idmapping_job(["P11111"], RUN_URL)

    def test_single_accession(self):
        resp = _mock_response({"jobId": "solo"})
        with patch(f"{MODULE}.requests.post", return_value=resp):
            job_id = _submit_idmapping_job(["P11111"], RUN_URL)
        assert job_id == "solo"


# ---------------------------------------------------------------------------
# _wait_for_results
# ---------------------------------------------------------------------------

class TestWaitForResults:
    def test_returns_immediately_when_results_present(self):
        resp = _mock_response({"results": [{"from": "P11111", "to": "UPI000"}]})
        with patch(f"{MODULE}.requests.get", return_value=resp), \
             patch(f"{MODULE}.time.sleep") as mock_sleep:
            _wait_for_results("job1", RESULTS_URL, poll_interval=0)
        mock_sleep.assert_not_called()

    def test_polls_until_results_appear(self):
        not_ready = _mock_response({"failedIds": []})   # no "results" key
        ready = _mock_response({"results": []})
        with patch(f"{MODULE}.requests.get", side_effect=[not_ready, not_ready, ready]), \
             patch(f"{MODULE}.time.sleep") as mock_sleep:
            _wait_for_results("job1", RESULTS_URL, poll_interval=0.01)
        assert mock_sleep.call_count == 2

    def test_uses_correct_url(self):
        resp = _mock_response({"results": []})
        with patch(f"{MODULE}.requests.get", return_value=resp) as mock_get, \
             patch(f"{MODULE}.time.sleep"):
            _wait_for_results("job42", RESULTS_URL, poll_interval=0)
        called_url = mock_get.call_args[0][0]
        assert "job42" in called_url

    def test_sleeps_with_given_interval(self):
        not_ready = _mock_response({})
        ready = _mock_response({"results": []})
        with patch(f"{MODULE}.requests.get", side_effect=[not_ready, ready]), \
             patch(f"{MODULE}.time.sleep") as mock_sleep:
            _wait_for_results("job1", RESULTS_URL, poll_interval=0.5)
        mock_sleep.assert_called_once_with(0.5)

    def test_non_200_response_keeps_polling(self):
        error_resp = _mock_response(status_code=503)
        error_resp.raise_for_status.return_value = None  # _wait_for_results doesn't call raise_for_status
        ready = _mock_response({"results": []})
        with patch(f"{MODULE}.requests.get", side_effect=[error_resp, ready]), \
             patch(f"{MODULE}.time.sleep"):
            _wait_for_results("job1", RESULTS_URL, poll_interval=0)  # should not raise


# ---------------------------------------------------------------------------
# _fetch_paginated_results
# ---------------------------------------------------------------------------

class TestFetchPaginatedResults:
    def test_returns_mapping_from_single_page(self):
        payload = {"results": [{"from": "P11111", "to": "UPI001"}, {"from": "P22222", "to": "UPI002"}]}
        resp = _mock_response(payload, links={})
        with patch(f"{MODULE}.requests.get", return_value=resp):
            result = _fetch_paginated_results("job1", RESULTS_URL)
        assert result == {"P11111": "UPI001", "P22222": "UPI002"}

    def test_follows_next_link_across_pages(self):
        page1 = _mock_response(
            {"results": [{"from": "P11111", "to": "UPI001"}]},
            links={"next": {"url": "https://next-page.url"}},
        )
        page2 = _mock_response(
            {"results": [{"from": "P22222", "to": "UPI002"}]},
            links={},
        )
        with patch(f"{MODULE}.requests.get", side_effect=[page1, page2]):
            result = _fetch_paginated_results("job1", RESULTS_URL)
        assert result == {"P11111": "UPI001", "P22222": "UPI002"}

    def test_stops_when_no_next_link(self):
        resp = _mock_response({"results": [{"from": "P11111", "to": "UPI001"}]}, links={})
        with patch(f"{MODULE}.requests.get", return_value=resp) as mock_get:
            _fetch_paginated_results("job1", RESULTS_URL)
        assert mock_get.call_count == 1

    def test_empty_results_returns_empty_dict(self):
        resp = _mock_response({"results": []}, links={})
        with patch(f"{MODULE}.requests.get", return_value=resp):
            result = _fetch_paginated_results("job1", RESULTS_URL)
        assert result == {}

    def test_raises_on_http_error(self):
        resp = _mock_response(status_code=500)
        with patch(f"{MODULE}.requests.get", return_value=resp):
            with pytest.raises(requests.HTTPError):
                _fetch_paginated_results("job1", RESULTS_URL)

    def test_merges_results_across_three_pages(self):
        def make_page(pairs, next_url=None):
            links = {"next": {"url": next_url}} if next_url else {}
            return _mock_response(
                {"results": [{"from": f, "to": t} for f, t in pairs]},
                links=links,
            )
        pages = [
            make_page([("P11111", "UPI001")], next_url="https://p2"),
            make_page([("P22222", "UPI002")], next_url="https://p3"),
            make_page([("P33333", "UPI003")]),
        ]
        with patch(f"{MODULE}.requests.get", side_effect=pages):
            result = _fetch_paginated_results("job1", RESULTS_URL)
        assert result == {"P11111": "UPI001", "P22222": "UPI002", "P33333": "UPI003"}

    def test_uses_job_id_in_initial_url(self):
        resp = _mock_response({"results": []}, links={})
        with patch(f"{MODULE}.requests.get", return_value=resp) as mock_get:
            _fetch_paginated_results("jobXYZ", RESULTS_URL)
        first_url = mock_get.call_args_list[0][0][0]
        assert "jobXYZ" in first_url
