import logging
from io import StringIO
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import requests
import pandas as pd
from clans3d.core.input_file_type import InputFileType
from clans3d.utils.api_utils import uniprot_accessions_to_uniparc_accessions

logger = logging.getLogger(__name__)


def extract_uids_from_fasta(fasta_file):
    """
    Extract all UniProt IDs from Fasta file.
    All entries must contain a Uniprot_ID in the header.
    :return: a list which contains all Uniprot_IDs
    """
    uids = []
    for entry in SeqIO.parse(fasta_file, "fasta"):
        header = entry.id
        uid = extract_uid_from_recordID(header)
        uids.append(uid)
    return uids


def clean_fasta_file(fasta_file, uids):
    """
    Cleans the fasta file by removing entries that are not in the given uid list.
    :param fasta_file: path to the input fasta file
    :param uids: list of UniProt IDs to keep
    :return: path to the cleaned fasta file
    """
    cleaned_fasta_path = "cleaned.fasta"
    with open(cleaned_fasta_path, 'w') as cleaned_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in uids:
                SeqIO.write(record, cleaned_fasta, "fasta")
    return cleaned_fasta_path


def extract_uid_from_recordID(record_id):
    """
    Extracts the UniProt ID from a record ID.
    It removes any version numbers and prefixes.
    Assumes the record ID is in the format 'tr|A0A2M8UWB6.1|A0A2M8UWB6_PSESP/524-595'.
    :param record_id: the record ID string
    :return: the extracted UniProt ID
    """
    if "|" in record_id:
        uid = record_id.split("|")[1].split(".")[0]
    elif "_" in record_id:
        uid = record_id.split("_")[0].split(".")[0]
    elif "/" in record_id:
        uid = record_id.split("/")[0].split(".")[0]
    else:
        # if no '|' is found in header, the header is assumed to be the uid
        uid = record_id.split(".")[0]
    return uid


def extract_region_from_record(record) -> tuple[int, int] | None:
    """
    Extracts the numeric region of interest (start and end) from a FASTA record header.

    Example header:
        >Y1502_ARCFU/1-68
        >tr|A0A819D2C1|A0A819D2C1_9BILA/1286-1334

    This function extracts the region "1-68" or "1286-1334" and returns [start, end] as integers.

    Args:
        record (SeqRecord): A Biopython SeqRecord object representing a FASTA entry.

    Returns:
        list[int, int] | None: A tuple (region_start, region_end) if a valid region is found,
                               otherwise None.
    """
    # Extract the record header (identifier)
    header = record.id
    # Search for a region pattern like "/123-456"
    match = re.search(r"/(\d+)-(\d+)", header)
    if match:
        start, end = map(int, match.groups())
        return (start, end)
    else:
        return None


def add_region_to_record(record: SeqRecord, region: tuple[int, int] | None) -> SeqRecord:
    """
    Takes a SeqRecord and modifies its ID and Sequence to include the region if provided.

    Args:
        record (SeqRecord): The record to modify.
        region (tuple[int, int] | None): The region to include in the ID and to cut the sequence with.

    Returns:
        SeqRecord: The modified record.
    """
    if region is None:
        region_str = ""
    else:
        if record.seq is None:
            raise ValueError(f"Record {record.id} has no sequence to extract region from.")
        seq_len = len(record.seq)
        start, end = region
        if start < 1 or end > seq_len or start > end:
            raise ValueError(
                f"Invalid region {region} for record {record.id} with length {seq_len}")
            
        region_str = f"/{int(region[0])}-{int(region[1])}"
        record = record[region[0]-1 : region[1]]
    record.id = f"{record.id}{region_str}"
    record.description = ""
    return record


def create_mock_up_record(uid: str, region: tuple[int, int] | None) -> SeqRecord:
    """
    Creates a mock-up SeqRecord for a given UID and region.

    Args:
        uid (str): The UniProt ID.
        region (tuple[int, int] | None): The region to include in the ID.

    Returns:
        SeqRecord: The created mock-up record.
    """
    if region is None:
        region_str = ""
    else:
        region_str = f"/{int(region[0])}-{int(region[1])}"
    record_id = f"|{uid}|{region_str}"
    record = SeqRecord(Seq("not_found"), id=record_id, description="")
    return record


def copy_records_from_fasta(path_to_fasta: str, uids: list[str], out_path: str) -> str:
    """
    Copies the records with the given uids from the fasta file to a new fasta file.
    Args:
        path_to_fasta (str): Path to the input FASTA file.
        uids (list[str]): List with UIDs to copy.
        out_path (str): Path to the output FASTA file.

    Returns:
        str: Path to the output FASTA file.
    """
    original_records = list(SeqIO.parse(path_to_fasta, "fasta"))
    records = []
    for uid in uids:
        for record in original_records:
            if uid in record.id:
                records.append(record)
    SeqIO.write(records, out_path, "fasta")
    return out_path


def download_fasta_record(uid: str, upi=None, region: tuple[int, int] | None = None) -> SeqRecord | bool:
    """
    Download a FASTA record by accession from Uniprot and as a fallback from UniParc.
    If the record is downloaded from UniParc the initial uid is used as accession-id.
    If the sequence is not found, it returns False.
    If a region is provided, the sequence is truncated to that region and the record ID is modified accordingly.   

    Args:
        uid (str): UniProt accession.
        upi (str, optional): UniParc accession. Defaults to None.
        region (tuple[int, int] | None): Residue range [start, end]
    Returns:
        SeqRecord: downloaded FASTA record.
        bool: False if download failed.
    """
    logger.debug("Attempting to download FASTA record for UID: %s with UPI: %s and region: %s", uid, upi, region)
    uniprot_api_url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
    uniparc_api_url = f"https://rest.uniprot.org/uniparc/{upi}.fasta"
    urls = [uniprot_api_url, uniparc_api_url if upi is not None else uniprot_api_url]
    for url in urls:
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200 and response.text.startswith(">"):
                fasta_io = StringIO(response.text)
                record = next(SeqIO.parse(fasta_io, "fasta"))
                if region:
                    record = add_region_to_record(record, region)
                return record
        except Exception as e:
            continue
    return False


def remove_non_existing_uniprot_accessions(uniprot_accessions: list[str]) -> list[str]:
    """
    Given a list of uniprot accessions, the method removes those which do not have an existing Uniprot or Uniparc entry.

    Args:
        uniprot_accessions (list[str]): List of uniprot accessions.

    Returns:
        list[str]: Filtered list of uniprot accessions with existing entries.
    """
    uids_to_upis = uniprot_accessions_to_uniparc_accessions(uniprot_accessions)
    for uid, upi in uids_to_upis.items():
        record = download_fasta_record(uid, upi=upi)
        if record is False:
            logger.warning("Could not retrieve sequence for UID: %s", uid)
            uniprot_accessions.remove(uid)
    return uniprot_accessions


def filter_input_file(in_path: str, out_path: str, dataset_type: InputFileType):
    """
    Given a path to a FASTA file or TSV file, the method filters out entries without existing Uniprot or Uniparc entries. 
    It saves the filtered dataset to out_path. The TSV file must contain a column 'entry' with the Uniprot accessions.

    Args:
        in_path (str): Path to the input dataset file.
        out_path (str): Path to the output file.
        dataset_type (InputFileType): Type of the input dataset file. Can be FASTA or TSV.
    """
    logger.debug("Filtering...")
    if dataset_type == InputFileType.FASTA:
        logger.debug("Extracting uids...")
        uids = extract_uids_from_fasta(in_path)
        uids_filtered = remove_non_existing_uniprot_accessions(uids)
        copy_records_from_fasta(in_path, uids_filtered, out_path)
    elif dataset_type == InputFileType.TSV:
        logger.debug("Extracting uids...")
        df = pd.read_csv(in_path, sep='\t')
        uids = df['entry'].tolist()
        uids_filtered = remove_non_existing_uniprot_accessions(uids)
        df_filtered = df[df['entry'].isin(uids_filtered)]
        df_filtered.to_csv(out_path, sep='\t', index=False)
    else:   
        raise ValueError("Unsupported dataset type. Use InputFileType.FASTA or InputFileType.TSV.")


def generate_fasta_from_uids_with_regions(uids_with_regions: dict[str, tuple[int, int] | None], out_path: str, original_fasta=None):
    """
    Generates a FASTA file containing the sequences (cut down to their corresponding regions) of the given UIDs.
    The records will include the UID and region (if present) in the header.
    
    If original_fasta is provided, the method will copy sequences which are also part of the uids_with_regions.
    Otherwise the fasta file is generated from scratch and the Record IDs are of the format: UID or UID/start-end.
    The sequences are cut to the specified regions if provided.

    Args:
        uids_with_regions (dict[str, tuple[int, int] | None]): Mapping {uid: (region_start, region_end) | None}.
        out_path (str): Path to the output FASTA file.
        original_fasta (str, optional): Path to an existing FASTA file with sequences.
    Returns:
    """
    logger.debug("Generating FASTA from UIDs with regions...")
    uids = list(uids_with_regions.keys())
    # generate fasta with records of original_fasta file
    if original_fasta is not None:
        copy_records_from_fasta(original_fasta, uids, out_path)
    # generate fasta with possible regions from scratch
    else:
        records = []
        uids_to_upis = uniprot_accessions_to_uniparc_accessions(uids)
        for uid, region in uids_with_regions.items():
            upi = uids_to_upis.get(uid)
            downloaded_record_with_region = download_fasta_record(uid, upi, region)
            if isinstance(downloaded_record_with_region, SeqRecord):
                downloaded_record_with_region.id = uid
                records.append(downloaded_record_with_region)
            else:
                fallback_record = create_mock_up_record(uid, region)
                records.append(fallback_record)
        SeqIO.write(records, out_path, "fasta")
    return out_path
