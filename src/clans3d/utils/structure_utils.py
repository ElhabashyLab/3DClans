import logging
import os
import requests
import pandas as pd
from Bio import SeqIO
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.mmcifio import MMCIFIO
from Bio.PDB.PDBIO import Select
from clans3d.core.input_file_type import InputFileType
from clans3d.utils.file_utils import download_file, reset_dir_content
from clans3d.utils.fasta_utils import extract_uid_from_recordID, extract_region_from_record

logger = logging.getLogger(__name__)


def _log_interval(total: int) -> int:
    """Return how often to log progress: every 10 items or every 20%, whichever is smaller."""
    return max(1, min(10, total // 5))


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


def _fetch_alphafold_cif_url(accession_id: str) -> str | None:
    """
    Queries the AlphaFold DB REST API and returns the CIF download URL for the
    given accession, or None if the request fails or no URL is found.

    Args:
        accession_id (str): UniProt accession to look up.

    Returns:
        str | None: CIF download URL, or None on failure.
    """
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{accession_id}"
    try:
        response = requests.get(api_url, timeout=30)
        response.raise_for_status()
        predictions = response.json()
    except requests.RequestException as e:
        logger.debug("Failed to reach AlphaFold API for %s: %s", accession_id, e)
        return None
    if not predictions:
        logger.debug("No predictions found in AlphaFold API response for %s.", accession_id)
        return None
    url = predictions[0].get("cifUrl")
    if not url:
        logger.debug("No cifUrl in AlphaFold API response for %s.", accession_id)
    return url


def fetch_structures(input_file_path: str, input_file_type: InputFileType, output_dir: str) -> dict[str, tuple[int, int] | None]:
    """
    Fetches and stores structure files, specified in a given input_file in a output directory.
    The contents of the output dir will be overwritten.
    It returns a dict containing {uid : (region_start, region_end)} of the sequences that have been successfully downloaded.
    
    Args:
        input_file_path: Path to the input file.
        input_file_type: Specifies the type of the input file.
        output_dir: The directory where the downloaded structures are stored.
        
    Returns:
        dict (str, tuple[int, int] | None): Containing downloaded uids together with their regions.
    """    
    reset_dir_content(output_dir)
    if input_file_type is InputFileType.FASTA or input_file_type is InputFileType.A2M:
        return process_fasta_file(input_file_path, output_dir)
    elif input_file_type is InputFileType.TSV:
        return process_tsv_file(input_file_path, output_dir)
    else:
        raise ValueError(f"Unsupported input file type: {input_file_type}")
    

def process_fasta_file(input_file_path: str, output_dir: str) -> dict:
    """
    Downloads the structures of the sequences in the given fasta_file.
    It also creates a dictionary with an entry for each downloaded record: {uid : region}
    The region can be None if not specified in the header of the fasta records.

    Args:
        input_file_path (str): Path to fasta file.
        output_dir (str): Directory in which to save the downloaded structures.

    Returns:
        dict: Containing downloaded uids together with their regions.
    """
    successful_downloads = 0
    uids_with_regions = {}
    records = list(SeqIO.parse(input_file_path, "fasta"))
    total_records = len(records)
    interval = _log_interval(total_records)
    for idx, record in enumerate(records, 1):
        try:
            uid = extract_uid_from_recordID(record.id)
            region = extract_region_from_record(record)
        except ValueError as e:
            logger.warning("Skipping entry '%s': %s", record.id, e)
            if idx % interval == 0 or idx == total_records:
                logger.info("Downloaded %d/%d structures...", idx, total_records)
            continue
        if download_alphafold_structure(uid, output_dir, region):
            uids_with_regions[uid] = region
            successful_downloads += 1
        else:
            logger.warning("Skipping entry '%s': could not download AlphaFold structure.", uid)
        if idx % interval == 0 or idx == total_records:
            logger.info("Downloaded %d/%d structures...", idx, total_records)
    failed = total_records - successful_downloads
    logger.info("Downloaded %d/%d structure files (%d failed).", successful_downloads, total_records, failed)
    return uids_with_regions


def process_tsv_file(input_file_path: str, output_dir: str) -> dict:
    """
    Downloads the structures of the uids in the given tsv file.
    It also creates a dictionary with an entry for each downloaded record. {uid : [region]}
    The tsv file should contain the columns [entry, region_start, region_end].

    Args:
        input_file_path (str): Path to tsv file.
        output_dir (str): Directory in which to save the downloaded structures.

    Returns:
        dict: Containing downloaded records together with their regions.
    """
    df = pd.read_csv(input_file_path, sep='\t')
    successful_downloads = 0
    uids_with_regions = {}
    total_uids = len(df)
    interval = _log_interval(total_uids)
    for idx, (_, row) in enumerate(df.iterrows(), 1):
        uid = row['entry']
        try:
            region = _parse_region_from_tsv_row(row)
        except ValueError as e:
            logger.warning("Skipping entry '%s': %s", uid, e)
            continue
        if download_alphafold_structure(uid, output_dir, region):
            successful_downloads += 1
            uids_with_regions[uid] = region
        else:
            logger.warning("Skipping entry '%s': could not download AlphaFold structure.", uid)
        if idx % interval == 0 or idx == total_uids:
            logger.info("Downloaded %d/%d structures...", idx, total_uids)
    failed = total_uids - successful_downloads
    logger.info("Downloaded %d/%d structure files (%d failed).", successful_downloads, total_uids, failed)
    return uids_with_regions


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
