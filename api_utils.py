import requests
import time
from typing import Iterable


def _chunked(iterable: list[str], size: int) -> Iterable[list[str]]:
    """
    Yield successive chunks from a list.    
    Args:
        iterable (List[str]): Input list.
        size (int): Maximum chunk size.

    Yields:
        List[str]: Chunks of the input list.
    """
    for i in range(0, len(iterable), size):
        yield iterable[i : i + size]


def _submit_idmapping_job(uniprot_accessions: list[str], UNIPROT_IDMAPPING_RUN_URL: str) -> str:
    """
    Submit a UniProt ID mapping job.

    Args:
        uniprot_accessions (list[str]): UniProt accessions to map.
        UNIPROT_IDMAPPING_RUN_URL (str): Url to uniprot uid to upi mapping.

    Returns:
        str: UniProt job ID.
    """
    response = requests.post(
        UNIPROT_IDMAPPING_RUN_URL,
        data={
            "from": "UniProtKB_AC-ID",
            "to": "UniParc",
            "ids": ",".join(uniprot_accessions),
        },
    )
    response.raise_for_status()
    return response.json()["jobId"]


def _wait_for_results(job_id: str, UNIPROT_IDMAPPING_RESULTS_URL: str, poll_interval: float = 1.0) -> None:
    """
    Block until UniProt mapping results are available.

    UniProt does not provide a stable status schema, so completion is detected
    by polling the results endpoint until a response containing a "results"
    field is returned.

    Args:
        job_id (str): UniProt job ID.
        UNIPROT_IDMAPPING_RESULTS_URL (str): URL to results of IDMAPPING
        poll_interval (float): Seconds between polling attempts.
    """
    results_url = UNIPROT_IDMAPPING_RESULTS_URL.format(job_id)
    while True:
        response = requests.get(results_url)
        if response.status_code == 200:
            payload = response.json()
            if "results" in payload:
                return
        time.sleep(poll_interval)


def _fetch_paginated_results(job_id: str, UNIPROT_IDMAPPING_RESULTS_URL: str) -> dict[str, str]:
    """
    Fetch all paginated mapping results for a completed job.

    Args:
        job_id (str): UniProt job ID.
        UNIPROT_IDMAPPING_RESULTS_URL (str): URL to results of IDMAPPING

    Returns:
        dict[str, str]: Mapping from UniProt accession to UniParc accession.
    """
    results_url = UNIPROT_IDMAPPING_RESULTS_URL.format(job_id)
    mapping: dict[str, str] = {}
    next_url = results_url

    while next_url:
        response = requests.get(next_url)
        response.raise_for_status()
        payload = response.json()

        for item in payload.get("results", []):
            mapping[item["from"]] = item["to"]

        next_url = response.links.get("next", {}).get("url")

    return mapping


def uniprot_accessions_to_uniparc_accessions(
    uniprot_accessions: list[str],
    batch_size: int = 200,
) -> dict[str, str | None]:
    """
    Convert UniProt accessions to UniParc accessions using the UniProt
    ID mapping service.

    Uses batching and returns None for accessions that were not mapped

    Args:
        uniprot_accessions (list[str]): UniProt accessions to convert.
        batch_size (int): Maximum number of accessions per UniProt job.

    Returns:
        dict[str, str | None]: Mapping {uniprot_accession: uniparc_accession | None}.
    """
    UNIPROT_IDMAPPING_RUN_URL = "https://rest.uniprot.org/idmapping/run"
    UNIPROT_IDMAPPING_RESULTS_URL = "https://rest.uniprot.org/idmapping/results/{}"
    final_mapping: dict[str, str | None] = {}
    for batch in _chunked(uniprot_accessions, batch_size):
        job_id = _submit_idmapping_job(batch, UNIPROT_IDMAPPING_RUN_URL)
        _wait_for_results(job_id, UNIPROT_IDMAPPING_RESULTS_URL)
        batch_mapping = _fetch_paginated_results(job_id, UNIPROT_IDMAPPING_RESULTS_URL)
        for acc in batch:
            final_mapping[acc] = batch_mapping.get(acc)
    return final_mapping


def uniprot_accessions_to_uniparc_accessions_old(uniprot_accessions: list[str]) -> dict:
    """
    Converts a UniProt accessions to UniParc accessions using the UniProt mapping service.

    Args:
        uniprot_accessions (lst): The UniProt accessions to convert.

    Returns:
        dict: containing the mapping {uniprot_accession: uniparc_accession|None}.
    """
    url = "https://rest.uniprot.org/idmapping/run"
    params = {
        "from": "UniProtKB_AC-ID",
        "to": "UniParc",
        "ids": ",".join(uniprot_accessions)
    }
    response = requests.post(url, data=params)
    response.raise_for_status()
    job_id = response.json()["jobId"]
    # Wait for completion
    result_url = f"https://rest.uniprot.org/idmapping/results/{job_id}"
    status_url = f"https://rest.uniprot.org/idmapping/status/{job_id}"
    while True:
        r = requests.get(result_url)
        if r.status_code == 200 and "results" in r.json():
            break
        time.sleep(1)
    # Fetch paginated results
    mapping = {}
    next_url = result_url
    while next_url:
        r = requests.get(next_url)
        r.raise_for_status()
        data = r.json()
        for item in data.get("results", []):
            mapping[item["from"]] = item["to"]
        next_url = r.links.get("next", {}).get("url")
    # Fill missing with None
    return {acc: mapping.get(acc) for acc in uniprot_accessions}
