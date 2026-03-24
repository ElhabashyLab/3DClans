import logging
import os
import requests
from typing import Literal
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
from Bio import SeqIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import Select
from clans3d.core.input_file_type import InputFileType
from clans3d.utils.file_utils import download_file, reset_dir_content
from clans3d.utils.fasta_utils import extract_uid_from_recordID, extract_region_from_record
from clans3d.utils.log import _log_interval

logger = logging.getLogger(__name__)


def _parse_region_from_tsv_row(row: pd.Series) -> tuple[int, int] | None:
    """
    Parses and validates the region columns from a TSV row.

    Returns:
        tuple[int, int]: (start, end) if both are present and valid.
        None: if both are absent (deliberately no region specified).

    Raises:
        ValueError: if only one of start/end is present, or if the range is biologically invalid.
    """
    start = row.get('region_start', pd.NA)
    end = row.get('region_end', pd.NA)
    start_missing = pd.isna(start)
    end_missing = pd.isna(end)
    if start_missing and end_missing:
        return None
    if start_missing or end_missing:
        raise ValueError("only one of region_start/region_end is specified")
    start, end = int(start), int(end)
    if start < 1 or end < 1 or end < start:
        raise ValueError(f"invalid region (start={start}, end={end})")
    return (start, end)


def _fetch_alphafold_cif_url(accession_id: str) -> str | Literal[False]:
    """
    Queries the AlphaFold DB REST API and returns the CIF download URL for the
    given accession, or False if the request fails or no URL is found.

    Args:
        accession_id (str): UniProt accession to look up.

    Returns:
        str | Literal[False]: CIF download URL, or False on failure.
    """
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession_id}"
    try:
        response = requests.get(api_url, timeout=30)
        response.raise_for_status()
        predictions = response.json()
    except requests.RequestException as e:
        logger.debug("Failed to reach AlphaFold API for %s: %s", accession_id, e)
        return False
    if not predictions:
        logger.debug("No predictions found in AlphaFold API response for %s.", accession_id)
        return False
    url = predictions[0].get("cifUrl")
    if not url:
        logger.debug("No cifUrl in AlphaFold API response for %s.", accession_id)
        return False
    return url


def _run_download_tasks(
    tasks: list[tuple[str, tuple[int, int] | None]],
    output_dir: str,
    max_workers: int = 10,
) -> dict[str, tuple[int, int] | None]:
    """Run download/optional-region-extraction tasks concurrently.

    Args:
        tasks: List of ``(uid, region)`` tuples to process.
        output_dir: Directory where structures are written.
        max_workers: Number of concurrent download threads.

    Returns:
        Mapping of successfully downloaded UIDs to their regions.
    """
    total = len(tasks)
    interval = _log_interval(total)
    successful_downloads = 0
    uids_with_regions: dict[str, tuple[int, int] | None] = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_task = {
            executor.submit(download_alphafold_structure, uid, output_dir, region): (uid, region)
            for uid, region in tasks
        }
        completed = 0
        for future in as_completed(future_to_task):
            uid, region = future_to_task[future]
            completed += 1
            if future.result():
                uids_with_regions[uid] = region
                successful_downloads += 1
            else:
                logger.warning("Skipping entry '%s': could not download AlphaFold structure.", uid)
            if completed % interval == 0 or completed == total:
                logger.info("Processing %d/%d structures...", completed, total)

    failed = total - successful_downloads
    logger.info("Downloaded %d/%d structure files (%d failed).", successful_downloads, total, failed)
    return uids_with_regions


def fetch_structures(
    input_file_path: str,
    input_file_type: InputFileType,
    output_dir: str,
    max_workers: int = 10,
) -> dict[str, tuple[int, int] | None]:
    """
    Fetches and stores structure files, specified in a given input_file in a output directory.
    The contents of the output dir will be overwritten.
    It returns a dict containing {uid : (region_start, region_end)} of the sequences that have been successfully downloaded.
    Threading is used to speed up the process of downloading and region extraction.
    ``max_workers`` controls the pool size
    
    Args:
        input_file_path: Path to the input file.
        input_file_type: Specifies the type of the input file.
        output_dir: The directory where the downloaded structures are stored.
        max_workers: Maximum number of concurrent download threads (default: 10).
        
    Returns:
        dict (str, tuple[int, int] | None): Containing downloaded uids together with their regions.
    """    
    reset_dir_content(output_dir)
    if input_file_type in (InputFileType.FASTA, InputFileType.A2M, InputFileType.A3M):
        return process_fasta_file(input_file_path, output_dir, max_workers=max_workers)
    elif input_file_type is InputFileType.TSV:
        return process_tsv_file(input_file_path, output_dir, max_workers=max_workers)
    else:
        raise ValueError(f"Unsupported input file type: {input_file_type}")
    

def process_fasta_file(input_file_path: str, output_dir: str, max_workers: int = 10) -> dict:
    """
    Downloads the structures of the sequences in the given fasta_file.
    It also creates a dictionary with an entry for each downloaded record: {uid : region}
    The region can be None if not specified in the header of the fasta records.

    Downloads and region extraction run concurrently (``max_workers`` threads).
    Results are collected in order of completion; progress is logged from
    the main thread so no locking is needed.

    Args:
        input_file_path (str): Path to fasta file.
        output_dir (str): Directory in which to save the downloaded structures.
        max_workers (int): Number of concurrent download threads (default: 10).

    Returns:
        dict: Containing downloaded uids together with their regions.
    """
    records = list(SeqIO.parse(input_file_path, "fasta"))

    # Parse all (uid, region) pairs first so we can report skips before submitting
    tasks: list[tuple[str, tuple[int, int] | None]] = []
    for record in records:
        try:
            uid = extract_uid_from_recordID(record.id)
            region = extract_region_from_record(record)
            tasks.append((uid, region))
        except ValueError as e:
            logger.warning("Skipping entry '%s': %s", record.id, e)

    return _run_download_tasks(tasks, output_dir, max_workers=max_workers)


def process_tsv_file(input_file_path: str, output_dir: str, max_workers: int = 10) -> dict:
    """
    Downloads the structures of the uids in the given tsv file.
    It also creates a dictionary with an entry for each downloaded record. {uid : [region]}
    The tsv file should contain the columns [entry, region_start, region_end].

    Downloads and region extraction run concurrently (``max_workers`` threads).
    Results are collected in order of completion; progress is logged from
    the main thread so no locking is needed.

    Args:
        input_file_path (str): Path to tsv file.
        output_dir (str): Directory in which to save the downloaded structures.
        max_workers (int): Number of concurrent download threads (default: 10).

    Returns:
        dict: Containing downloaded records together with their regions.
    """
    df = pd.read_csv(input_file_path, sep='\t')

    # Parse all (uid, region) pairs first so we can report skips before submitting
    tasks: list[tuple[str, tuple[int, int] | None]] = []
    for _, row in df.iterrows():
        uid = row['entry']
        try:
            region = _parse_region_from_tsv_row(row)
            tasks.append((uid, region))
        except ValueError as e:
            logger.warning("Skipping entry '%s': %s", uid, e)

    return _run_download_tasks(tasks, output_dir, max_workers=max_workers)


def download_alphafold_structure(
    accession_id: str,
    output_dir: str,
    region: tuple[int, int] | None,
) -> bool:
    """
    Downloads the AlphaFold-predicted CIF structure for a given accession and saves
    it to the output directory. If a region is provided, the structure is trimmed to
    that residue range after download.

    Args:
        accession_id (str): UniProt accession of the structure to download.
        output_dir (str): Directory where the structure file will be saved.
        region (tuple[int, int] | None): Residue range (start, end) to extract
            (1-based, inclusive). Pass None to keep the full structure.

    Returns:
        bool: True if the structure was successfully downloaded (and trimmed),
              False otherwise.
    """
    logger.debug("Downloading structure for %s with region: %s", accession_id, region)
    cif_url = _fetch_alphafold_cif_url(accession_id)
    if not cif_url:
        return False
    save_path = os.path.join(output_dir, f"{accession_id}.cif")
    if not download_file(cif_url, save_path):
        return False
    if region:
        extract_region_of_protein(save_path, region, save_path)
    return True
    

def extract_region_of_protein(path_to_protein: str, region: tuple[int, int], output_path: str | None = None) -> str:
    """
    Extracts a specific residue range from a CIF structure file and writes the
    filtered structure to a file.

    Args:
        path_to_protein (str): Path to the input CIF structure file.
        region (tuple[int, int]): (start, end) residue indices to extract (1-based, inclusive).
        output_path (str, optional): Destination path. Defaults to the input path with a '_roi' suffix.

    Returns:
        str: Path to the written output file.
    """
    logger.debug("Extracting region %s from protein file: %s", region, path_to_protein)
    class SelectRegion(Select):
        def accept_residue(self, residue):
            return region[0] <= residue.id[1] <= region[1]

    parser = MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", path_to_protein)
    if output_path is None:
        root, ext = os.path.splitext(path_to_protein)
        output_path = f"{root}_roi{ext}"
    io = MMCIFIO()
    io.set_structure(structure)
    io.save(output_path, select=SelectRegion())
    return output_path
