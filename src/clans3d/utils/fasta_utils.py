import logging
from concurrent.futures import ThreadPoolExecutor, as_completed
from io import StringIO
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import requests
import pandas as pd
from clans3d.core.input_file_type import InputFileType
from clans3d.utils.api_utils import uniprot_accessions_to_uniparc_accessions
from clans3d.utils.log import _log_interval

logger = logging.getLogger(__name__)

# Matches the official UniProt accession format:
# https://www.uniprot.org/help/accession_numbers
_UNIPROT_ACCESSION_RE = re.compile(
    r"[OPQ][0-9][A-Z0-9]{3}[0-9]|[A-NR-Z][0-9](?:[A-Z][A-Z0-9]{2}[0-9]){1,2}"
)

def extract_records_from_fasta(fasta_file):
    """
    Extract all records from Fasta file.
    All entries must contain a Uniprot_ID in the header.
    :return: a list which contains all records
    """
    records = []
    for entry in SeqIO.parse(fasta_file, "fasta"):
        header = entry.id
        uid = extract_uid_from_recordID(header)
        record = SeqRecord(entry.seq, id=uid, description="")
        records.append(record)
    return records


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


def clean_fasta_file(fasta_file: str, uids: list[str], out_path: str) -> str:
    """
    Cleans the fasta file by removing entries that are not in the given uid list
    and writes the result to *out_path*.

    Args:
        fasta_file (str): Path to the input FASTA file.
        uids (list[str]): UniProt IDs to retain. Matching is done against
            ``record.id`` as parsed by BioPython (the first whitespace-delimited
            token of the FASTA header).
        out_path (str): Destination path for the filtered FASTA file.

    Returns:
        str: *out_path* (the path passed in).
    """
    with open(out_path, 'w') as cleaned_fasta:
        for record in SeqIO.parse(fasta_file, "fasta"):
            if record.id in uids:
                SeqIO.write(record, cleaned_fasta, "fasta")
    return out_path


def extract_uid_from_recordID(record_id):
    """
    Extracts the UniProt accession from a record ID by searching for the
    standard UniProt accession pattern directly in the string.

    Handles formats such as:
        sp|P49811|MYOD1_PIG, tr|A0A2M8UWB6.1|A0A2M8UWB6_PSESP/524-595,
        P49811, P49811.1, P49811/1-300

    Falls back to stripping any version suffix (e.g. ``P49811.1`` → ``P49811``)
    or splitting on ``_`` (gene name style, e.g. ``MYOD1_PIG`` → ``MYOD1``) when
    no accession pattern is found.

    Args:
        record_id (str): the record ID string
    Returns:
        str: the extracted UniProt ID
    """
    match = _UNIPROT_ACCESSION_RE.search(record_id)
    if match:
        return match.group()
    else:
        raise ValueError(f"Could not extract UniProt accession from record ID: {record_id}")


def extract_region_from_record(record) -> tuple[int, int] | None:
    """
    Extracts the numeric region of interest (start and end) from a FASTA record header.
    If the region is invalid, it raises a ValueError.

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
        if start < 1 or end < 1 or end < start:
            raise ValueError(f"Invalid region {start}-{end} in record ID: {header}")
        region = (start, end)
    else:
        region = None
    return region


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
        region_str = f"{int(region[0])}-{int(region[1])}"
    record_id = f"{uid}/{region_str}"
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


def download_fasta_record(uid: str, upi: str | None = None, region: tuple[int, int] | None = None) -> SeqRecord:
    """
    Download a FASTA record by accession from UniProt, falling back to UniParc if not found.
    A region can optionally be specified to extract a subsequence from the full sequence.

    Args:
        uid (str): UniProt accession.
        upi (str | None): UniParc accession used as fallback. Defaults to None.
        region (tuple[int, int] | None): Residue range [start, end].
    Returns:
        SeqRecord: The downloaded FASTA record.
    Raises:
        requests.exceptions.RequestException: On network or transport failure.
        ValueError: If the sequence is not found on either endpoint.
    """
    logger.debug("Downloading FASTA record for UID: %s with UPI: %s and region: %s", uid, upi, region)

    response = requests.get(f"https://rest.uniprot.org/uniprotkb/{uid}.fasta", timeout=10)
    if not (response.status_code == 200 and response.text.startswith(">")):
        if upi is None:
            raise ValueError(f"Sequence not found for UID: {uid}")
        logger.debug("UID %s not found in UniProt, trying UniParc (%s)", uid, upi)
        response = requests.get(f"https://rest.uniprot.org/uniparc/{upi}.fasta", timeout=10)
        if not (response.status_code == 200 and response.text.startswith(">")):
            raise ValueError(f"Sequence not found for UID: {uid} (UniParc: {upi})")

    record = next(SeqIO.parse(StringIO(response.text), "fasta"))
    if region:
        record = add_region_to_record(record, region)
    return record


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
        try:
            download_fasta_record(uid, upi=upi)
        except (ValueError, requests.exceptions.RequestException) as e:
            logger.warning("Could not retrieve sequence for UID %s: %s — skipping", uid, e)
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



def _download_or_mock_fasta_record(
    uid: str,
    upi: str | None,
    region: tuple[int, int] | None,
) -> tuple[SeqRecord, bool]:
    """Download one FASTA record and fall back to a mock-up record on failure."""
    try:
        record = download_fasta_record(uid, upi, region)
        # Keep existing UID
        record.id = uid
        return record, True
    except (ValueError, requests.exceptions.RequestException) as e:
        logger.warning("Could not download %s: %s — using mock-up", uid, e)
        return create_mock_up_record(uid, region), False


def generate_fasta_from_uids_with_regions(
    uids_with_regions: dict[str, tuple[int, int] | None],
    out_path: str,
    max_workers: int = 10,
):
    """
    Generates a FASTA file by downloading sequences for the given UIDs and regions from UniProt.

    Downloads are performed concurrently with bounded workers while preserving
    deterministic output ordering based on the input mapping order.

    Args:
        uids_with_regions (dict[str, tuple[int, int] | None]): Mapping {uid: (region_start, region_end) | None}.
        out_path (str): Path to the output FASTA file.
        max_workers (int): Maximum number of concurrent download threads.

    Returns:
        str: Path to the output FASTA file.
    """

    logger.debug("Generating FASTA from UIDs with regions...")
    items = list(uids_with_regions.items())
    total = len(items)

    uids = [uid for uid, _ in items]
    uids_to_upis = uniprot_accessions_to_uniparc_accessions(uids)

    records: list[SeqRecord | None] = [None] * total
    successful_downloads = 0
    interval = _log_interval(total)

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        future_to_index = {
            executor.submit(
                _download_or_mock_fasta_record,
                uid,
                uids_to_upis.get(uid),
                region,
            ): idx
            for idx, (uid, region) in enumerate(items)
        }

        completed = 0
        for future in as_completed(future_to_index):
            idx = future_to_index[future]
            record, was_downloaded = future.result()
            records[idx] = record
            if was_downloaded:
                successful_downloads += 1

            completed += 1
            if completed % interval == 0 or completed == total:
                logger.info("Downloading FASTA records %d/%d...", completed, total)

    failed = total - successful_downloads
    logger.info(
        "Downloaded %d/%d FASTA records (%d mock-up).",
        successful_downloads,
        total,
        failed,
    )

    SeqIO.write([record for record in records if record is not None], out_path, "fasta")
    return out_path


def clean_aligned_sequence(seq: str) -> str:
    """
    Remove A2M/A3M gap characters from a sequence, keeping all amino acid residues.
    
    A2M/A3M format uses:
    - Uppercase letters: match states (residues aligned to query) - KEPT
    - Lowercase letters: insert states (insertions relative to query, but still real residues!) - KEPT
    - Dots (.): gaps in insert states - REMOVED
    - Dashes (-): gaps in match states - REMOVED
    
    This produces the TRUE sequence as it would appear in UniProt, not just the
    aligned portion. Both uppercase and lowercase letters are real amino acids.
    
    Args:
        seq: Raw sequence string from A2M/A3M file.
        
    Returns:
        Clean uppercase sequence with only gap characters removed.
    """
    return ''.join(c.upper() for c in seq if c.isalpha())


def generate_fasta_from_alignment_file(
    uids_with_regions: dict[str, tuple[int, int] | None],
    out_path: str,
    alignment_file: str
) -> str:
    """
    Generate clean FASTA from an A2M/A3M alignment file by stripping alignment characters.
    
    This function extracts sequences from the alignment file, removes lowercase insertions
    and gap characters, and writes clean uppercase sequences to the output file.
    Only sequences whose UIDs are in uids_with_regions (i.e., had successful structure downloads)
    are included.
    
    Args:
        uids_with_regions: Mapping {uid: (region_start, region_end) | None} of UIDs to include.
        out_path: Path to the output FASTA file.
        alignment_file: Path to the input A2M/A3M alignment file.
        
    Returns:
        Path to the output FASTA file.
    """
    logger.debug("Generating FASTA from alignment file with cleaned sequences...")
    records = []
    
    for record in SeqIO.parse(alignment_file, "fasta"):
        try:
            uid = extract_uid_from_recordID(record.id)
        except ValueError:
            continue  # Skip records with unextractable UIDs (these records were never downoaded)
            
        if uid not in uids_with_regions:
            continue  # Skip UIDs that weren't successfully downloaded
            
        clean_seq = clean_aligned_sequence(str(record.seq))
        region = uids_with_regions[uid]
        
        if region:
            record_id = f"{uid}/{region[0]}-{region[1]}"
        else:
            record_id = uid
            
        records.append(SeqRecord(Seq(clean_seq), id=record_id, description=""))
    
    SeqIO.write(records, out_path, "fasta")
    logger.debug("Wrote %d cleaned sequences to %s", len(records), out_path)
    return out_path
