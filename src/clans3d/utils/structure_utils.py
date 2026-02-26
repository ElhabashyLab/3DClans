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
    """Return how often to log progress: every 10 items or every 10%, whichever is smaller."""
    return max(1, min(10, total // 10))


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
        uid = extract_uid_from_recordID(record.id)
        region = extract_region_from_record(record)
        if download_alphafold_structure(uid, output_dir, region):
            uids_with_regions[uid] = region
            successful_downloads += 1
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
        start: int = row['region_start']
        end: int = row['region_end']
        if pd.isna(start) or pd.isna(end):
            region = None
        else:
            region = (int(start), int(end))
        if download_alphafold_structure(uid, output_dir, region):
            successful_downloads += 1
            uids_with_regions[uid] = region
        if idx % interval == 0 or idx == total_uids:
            logger.info("Downloaded %d/%d structures...", idx, total_uids)
    failed = total_uids - successful_downloads
    logger.info("Downloaded %d/%d structure files (%d failed).", successful_downloads, total_uids, failed)
    return uids_with_regions


def download_alphafold_structure(
    id: str,
    output_dir: str,
    region: tuple[int, int] | None,
    file_format: str = "cif",
) -> bool:
    """
    Downloads the AlphaFold-predicted structure for a given structure-ID (f.e. Uniprot-ID) and saves it
    to the specified output directory. If a residue region is provided, the structure
    is truncated to that region after download.
    The structure URL is obtained via the official AlphaFold DB REST API.

    Args:
        id (str): accession-id of the structure.
        output_dir (str): Directory where the structure file will be saved.
        region (tuple[int, int] | None): Residue range [start, end] to extract (1-based, inclusive).
        file_format (str): Structure format to download ("pdb" or "cif"). Defaults to "cif".

    Returns:
        bool: True if the structure was successfully downloaded (and processed),
              False otherwise.
    """
    logger.debug("Downloading structure for ID: %s with region: %s", id, region)
    if file_format not in {"pdb", "cif"}:
        raise ValueError("file_format must be 'pdb' or 'cif'")
    
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{id}"
    try:
        response = requests.get(api_url, timeout=30)
        response.raise_for_status()
        predictions = response.json()
    except Exception:
        logger.warning("Failed to download %s from %s", id, api_url)
        return False
    if not predictions:
        logger.warning("Failed to download %s from %s", id, api_url)
        return False
    
    model_info = predictions[0]
    if file_format == "pdb":
        structure_url = model_info.get("pdbUrl")
    else:
        structure_url = model_info.get("cifUrl")
    if not structure_url:
        logger.warning("Failed to download %s from %s", id, api_url)
        return False

    path_to_structure = os.path.join(output_dir, f"{id}.{file_format}")
    success = download_file(structure_url, path_to_structure)
    if not success:
        logger.warning("Failed to download %s from %s", id, api_url)
        return False

    if region:
        extract_region_of_protein(
            path_to_structure,
            file_format,
            region,
            path_to_structure,
        )
    return True
    

def extract_region_of_protein(path_to_protein: str, file_type: str, region: tuple[int, int], output_path = None):
    """
    Extracts a specific residue range from a protein structure file (PDB or CIF)
    and writes the filtered structure to a new file.
    
    Args:
        path_to_protein (str): Path to the protein structure file.
        file_type (str): 'pdb' or 'cif'.
        region (tuple[int, int]): (start, end) residue indices to extract (1-based, inclusive).
        output_path (str, optional): Path to save extracted structure. Defaults to adding '_roi' suffix.
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
