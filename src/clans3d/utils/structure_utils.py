import logging
import os
from shutil import copy2
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


def _build_uid_region_tasks(
    input_file_path: str,
    input_file_type: InputFileType,
) -> list[tuple[str, tuple[int, int] | None]]:
    """
    Parse input records into a list of (uid, region) tasks.

    Extracts UniProt accessions and optional region annotations from the input file.
    Logs warnings for any entries that fail to parse but continues processing.

    Args:
        input_file_path: Path to the input file (FASTA, A2M, A3M, or TSV).
        input_file_type: Type of the input file.

    Returns:
        list of (uid, region) tuples. Region is None if not specified in the input.

    Raises:
        ValueError: If the input file type is not supported.
    """
    if input_file_type in (InputFileType.FASTA, InputFileType.A2M, InputFileType.A3M):
        records = list(SeqIO.parse(input_file_path, "fasta"))
        tasks: list[tuple[str, tuple[int, int] | None]] = []
        for record in records:
            try:
                uid = extract_uid_from_recordID(record.id)
                region = extract_region_from_record(record)
                tasks.append((uid, region))
            except ValueError as e:
                logger.warning("Skipping entry '%s': %s", record.id, e)
        return tasks

    if input_file_type is InputFileType.TSV:
        df = pd.read_csv(input_file_path, sep='\t')
        tasks = []
        for _, row in df.iterrows():
            uid = row['entry']
            try:
                region = _parse_region_from_tsv_row(row)
                tasks.append((uid, region))
            except ValueError as e:
                logger.warning("Skipping entry '%s': %s", uid, e)
        return tasks

    raise ValueError(f"Unsupported input file type: {input_file_type}")


def _prepare_local_structure(
    uid: str,
    region: tuple[int, int] | None,
    structures_db: str,
    output_dir: str,
) -> bool:
    """
    Copy a local CIF structure to the working directory and optionally trim it by region.

    Copies a CIF file from the structures_db directory to the output directory.
    If a region is specified, extracts only the requested residue range from the
    copied structure file.

    Args:
        uid: UniProt accession (used to locate <uid>.cif in structures_db).
        region: Residue range (start, end) to extract, 1-based inclusive. None keeps full structure.
        structures_db: Directory containing source CIF files (<uid>.cif).
        output_dir: Directory where the prepared structure is written.

    Returns:
        True if the structure was successfully prepared; False if copy, extraction, or cleanup failed.
    """
    source_path = os.path.join(structures_db, f"{uid}.cif")
    target_path = os.path.join(output_dir, f"{uid}.cif")

    if not os.path.exists(source_path):
        logger.warning("Skipping entry '%s': local CIF not found at %s", uid, source_path)
        return False

    try:
        copy2(source_path, target_path)
    except OSError as e:
        logger.warning("Skipping entry '%s': could not copy local CIF from %s: %s", uid, source_path, e)
        return False

    if region is not None:
        try:
            extract_region_of_protein(target_path, region, target_path)
        except Exception as e:
            logger.warning("Skipping entry '%s': could not extract region %s from local CIF: %s", uid, region, e)
            try:
                os.remove(target_path)
            except OSError:
                pass
            return False

    return True


def _collect_prepared_structures(
    futures_to_indices: dict,
    tasks: list[tuple[str, tuple[int, int] | None]],
) -> dict[str, tuple[int, int] | None]:
    """
    Collect results from parallel local structure preparation tasks.

    Processes futures from a ThreadPoolExecutor submission dict, maps results back to
    (uid, region) pairs, tracks success/failure, and logs a summary.

    Args:
        futures_to_indices: Mapping of Future objects to their index in the task list.
        tasks: Original task list of (uid, region) tuples.

    Returns:
        Dictionary mapping successfully prepared UIDs to their regions.
    """
    total = len(tasks)
    interval = _log_interval(total)
    results: list[bool | None] = [None] * total
    successful = 0
    completed = 0

    for future in as_completed(futures_to_indices):
        idx = futures_to_indices[future]
        results[idx] = future.result()
        if results[idx]:
            successful += 1
        completed += 1
        if completed % interval == 0 or completed == total:
            logger.info("Preparing local structures %d/%d...", completed, total)

    uids_with_regions: dict[str, tuple[int, int] | None] = {}
    for idx, (uid, region) in enumerate(tasks):
        if results[idx]:
            uids_with_regions[uid] = region

    failed = total - successful
    logger.info(
        "Prepared %d/%d local structure files (%d failed).",
        successful,
        total,
        failed,
    )
    return uids_with_regions


def prepare_structures_from_local_db(
    input_file_path: str,
    input_file_type: InputFileType,
    structures_db: str,
    output_dir: str,
    max_workers: int = 10,
) -> dict[str, tuple[int, int] | None]:
    """
    Prepare local CIF structures for the pipeline.

    Parses the input file to extract UniProt accessions and optional region annotations.
    Copies matching CIF files from ``structures_db`` into ``output_dir`` and optionally
    trims them to the requested region. Preparation runs in parallel using a thread pool.

    The input file remains the source of truth for UIDs and optional regions.
    Matching CIF files must be named ``<uid>.cif`` in the structures_db directory.

    Args:
        input_file_path: Path to input file (FASTA, A2M, A3M, or TSV).
        input_file_type: Type of the input file.
        structures_db: Directory containing source CIF files (<uid>.cif).
        output_dir: Directory where prepared structures are written (will be reset).
        max_workers: Number of parallel preparation threads (default: 10).

    Returns:
        Dictionary mapping successfully prepared UIDs to their regions.
        UIDs for which the CIF was not found or region extraction failed are excluded.
    """
    reset_dir_content(output_dir)
    tasks = _build_uid_region_tasks(input_file_path, input_file_type)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_index = {
            executor.submit(_prepare_local_structure, uid, region, structures_db, output_dir): idx
            for idx, (uid, region) in enumerate(tasks)
        }
        uids_with_regions = _collect_prepared_structures(future_to_index, tasks)

    return uids_with_regions


def fetch_structures(
    input_file_path: str,
    input_file_type: InputFileType,
    output_dir: str,
    max_workers: int = 10,
) -> dict[str, tuple[int, int] | None]:
    """
    Download AlphaFold structures for sequences in an input file.

    Parses the input file (FASTA, A2M, A3M, or TSV) to extract UniProt accessions and
    optional region annotations, then downloads matching AlphaFold CIF structures and
    optionally trims them by region. The output directory is reset before download.
    Downloads and region extraction run in parallel using a thread pool.

    Args:
        input_file_path: Path to the input file.
        input_file_type: Type of the input file (FASTA, A2M, A3M, or TSV).
        output_dir: Directory where downloaded structures are stored (will be reset).
        max_workers: Number of concurrent download threads (default: 10).

    Returns:
        Dictionary mapping successfully downloaded UIDs to their regions.
        Region values are None if not specified in the input file.

    Raises:
        ValueError: If the input file type is unsupported.
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
    Download AlphaFold structures for sequences in a FASTA file.

    Parses the FASTA file to extract UniProt accessions and optional region annotations,
    then downloads matching AlphaFold CIF structures and optionally trims them by region.
    Downloads and region extraction run concurrently using a thread pool.

    Args:
        input_file_path: Path to FASTA file.
        output_dir: Directory in which to save the downloaded structures.
        max_workers: Number of concurrent download threads (default: 10).

    Returns:
        Dictionary mapping successfully downloaded UIDs to their regions.
        Region values are None if not specified in the FASTA header.
    """
    tasks = _build_uid_region_tasks(input_file_path, InputFileType.FASTA)
    return _run_download_tasks(tasks, output_dir, max_workers=max_workers)


def process_tsv_file(input_file_path: str, output_dir: str, max_workers: int = 10) -> dict:
    """
    Download AlphaFold structures for UIDs in a TSV file.

    Parses the TSV file to extract UniProt accessions and optional region annotations
    (columns: entry, region_start, region_end), then downloads matching AlphaFold CIF
    structures and optionally trims them by region. Downloads and region extraction
    run concurrently using a thread pool.

    Args:
        input_file_path: Path to TSV file with columns [entry, region_start, region_end].
        output_dir: Directory in which to save the downloaded structures.
        max_workers: Number of concurrent download threads (default: 10).

    Returns:
        Dictionary mapping successfully downloaded UIDs to their regions.
        Region values are None if not specified in the TSV file.
    """
    tasks = _build_uid_region_tasks(input_file_path, InputFileType.TSV)
    return _run_download_tasks(tasks, output_dir, max_workers=max_workers)


def download_alphafold_structure(
    accession_id: str,
    output_dir: str,
    region: tuple[int, int] | None,
) -> bool:
    """
    Download an AlphaFold-predicted CIF structure and optionally extract a region.

    Queries the AlphaFold DB REST API for the given accession, downloads the CIF file,
    and optionally trims it to the specified residue range.

    Args:
        accession_id: UniProt accession of the structure to download.
        output_dir: Directory where the structure file will be saved.
        region: Residue range (start, end) to extract, 1-based inclusive.
            Pass None to keep the full structure.

    Returns:
        True if the structure was successfully downloaded (and optionally trimmed),
        False if the download or region extraction failed.
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
    Extract a residue range from a CIF structure file.

    Parses a CIF structure, filters residues to the specified range, and writes
    the trimmed structure to an output file.

    Args:
        path_to_protein: Path to the input CIF structure file.
        region: Residue range (start, end) to extract, 1-based inclusive.
        output_path: Destination path. If None, appends '_roi' suffix to the input path.

    Returns:
        Path to the written output file.
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
