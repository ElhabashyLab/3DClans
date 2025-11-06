from io import StringIO
import re
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import Select
import requests
import os
import shutil
from InputFileType import InputFileType
import pandas as pd


def download_file(url, output_path):
    """
    Downloads a file from a given URL and saves it to a specified local path.
    """
    try:
        response = requests.get(url, stream=True, timeout=10)
        # raise exception if request failed
        response.raise_for_status()
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'wb') as f:
            # download file in 8KB chunks
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        return True
    except Exception as e:
        print(f"Failed to download {url}: {str(e)}")
        return False
    

def delete_dir_content(dir_path):
    """
    Deletes the content of the specified directory and creates it if it does not exists.
    """
    if os.path.exists(dir_path):
        shutil.rmtree(dir_path)
    os.makedirs(dir_path)


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
    else:
        # if no '|' is found in header, the header is assumed to be the uid
        uid = record_id.split(".")[0]
    return uid


def extract_region_from_record(record) -> list | None:
    """
    Extracts the numeric region of interest (start and end) from a FASTA record header.

    Example header:
        >Y1502_ARCFU/1-68
        >tr|A0A819D2C1|A0A819D2C1_9BILA/1286-1334

    This function extracts the region "1-68" or "1286-1334" and returns [start, end] as integers.

    Args:
        record (SeqRecord): A Biopython SeqRecord object representing a FASTA entry.

    Returns:
        list[int, int] | None: A list [region_start, region_end] if a valid region is found,
                               otherwise None.
    """
    # Extract the record header (identifier)
    header = record.id
    # Search for a region pattern like "/123-456"
    match = re.search(r"/(\d+)-(\d+)", header)
    if match:
        start, end = map(int, match.groups())
        return [start, end]
    else:
        return None


def fetch_pdbs(input_file_path: str, input_file_type: InputFileType, output_dir: str) -> dict:
    """
    Fetches and stores PDB files, specified in a given input_file in a output directory.
    The contents of the output dir will be overwritten.
    It returns a dict containing {uid : [region_start, region_end]} of the sequences that have been successfully downloaded.
    
    Args:
        input_file_path: Path to the input file.
        input_file_type: Specifies the type of the input file.
        output_dir: The directory where the downloaded PDBs are stored.
        
    Returns:
        dict: Path to cleaned version of the input_file containing only successfully downloaded elements.
    """    
    delete_dir_content(output_dir)
    if input_file_type is InputFileType.FASTA:
        return process_fasta_file(input_file_path, output_dir)
    elif input_file_type is InputFileType.TSV:
        return process_tsv_file(input_file_path, output_dir)
    else:
        raise ValueError(f"Unsupported input file type: {input_file_type}")
    

def process_fasta_file(input_file_path: str, output_dir: str) -> dict:
    """
    Downloads the PDBs of the sequences in the given fasta_file.
    It also creates a dictionary with an entry for each downloaded record: {uid : region}
    The region can be None if not specified in the header of the fasta records.

    Args:
        input_file_path (str): Path to fasta file.
        output_dir (str): Directory in which to save the downloaded PDBs.

    Returns:
        dict: Containing downloaded uids together with their regions.
    """
    successful_downloads = 0
    total_records = 0
    uids_with_regions = {}
    for record in SeqIO.parse(input_file_path, "fasta"):
        total_records += 1
        uid = extract_uid_from_recordID(record.id)
        region = extract_region_from_record(record)
        if download_alphafold_structure(uid, output_dir, region):
            uids_with_regions[uid] = region
            successful_downloads += 1
    print(f"Downloaded {successful_downloads} from {total_records} PDB files successfully.")
    return uids_with_regions


def process_tsv_file(input_file_path: str, output_dir: str) -> dict:
    """
    Downloads the PDBs of the uids in the given tsv file.
    It also creates a dictionary with an entry for each downloaded record. uid : [region]
    The tsv file should contain the columns [entry, region_start, region_end].

    Args:
        input_file_path (str): Path to tsv file.
        output_dir (str): Directory in which to save the downloaded PDBs.

    Returns:
        dict: Containing downloaded records together with their regions.
    """
    df = pd.read_csv(input_file_path, sep='\t')
    successful_downloads = 0
    uids_with_regions = {}
    total_uids = len(df)
    for _, row in df.iterrows():
        uid = row['entry']
        start: int = row['region_start']
        end: int = row['region_end']
        if pd.isna(start) or pd.isna(end):
            region = None
        else:
            region = [start, end]
        if download_alphafold_structure(uid, output_dir, region):
            successful_downloads += 1
            uids_with_regions[uid] = region
    print(f"Downloaded {successful_downloads} from {total_uids} PDB files successfully.")
    return uids_with_regions


def download_alphafold_structure(uid: str, output_dir: str, region: list[int] | None, file_format: str = "pdb") -> bool:
    """
    Downloads the AlphaFold structure of a given protein and saves it in the output_dir.
    If region is provided, the protein is cut to the specified region.

    Args:
        uid (str): Protein UniProt ID.
        output_dir (str): Directory to save the file.
        region (list[int] | None): Residue range [start, end].
        file_format (str): 'pdb' or 'cif'.
    Returns:
        bool: True if download succeeded.
    """
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v6.{file_format}" # This url might change over time depending on model
    path_to_structure = os.path.join(output_dir, f"{uid}.{file_format}")
    success = download_file(url, path_to_structure)
    if success and region:
        extract_region_of_protein(path_to_structure, file_format, region, path_to_structure) # overwrites path_to_structure
    return success
    

def extract_region_of_protein(path_to_protein: str, file_type: str, region: list[int], output_path = None):
    """
    Extracts a specific residue range from a protein structure file (PDB or CIF)
    without corrupting metadata. Writes output to a new file by default.
    
    Args:
        path_to_protein (str): Path to the protein file.
        file_type (str): 'pdb' or 'cif'.
        region (list[int]): [start, end] residue indices to extract.
        output_path (str, optional): Path to save extracted structure. Defaults to adding '_region' suffix.
    """
    class SelectRegion(Select):
        def accept_residue(self, residue):
            return region[0] <= residue.id[1] <= region[1]

    parser = PDBParser(QUIET=True) if file_type == "pdb" else MMCIFParser(QUIET=True)
    structure = parser.get_structure("protein", path_to_protein)
    if output_path is None:
        root, ext = os.path.splitext(path_to_protein)
        output_path = f"{root}_roi{ext}"
    meta_data = get_meta_data_of_structure_file(path_to_protein)
    meta_data_corrected = adapt_metadata_to_region(meta_data, region)
    with open(output_path, "w") as out:
        # out.writelines(meta_data_corrected)
        io = PDBIO()
        io.set_structure(structure)
        io.save(out, select=SelectRegion())
    return output_path


def adapt_metadata_to_region(meta_data: list[str], region: list[int]) -> list[str]:
    """
    Adapts the metadata lines of a structure file to reflect the extracted region.
    
    Args:
        meta_data (list[str]): Original metadata lines.
        region (list[int]): [start, end] residue indices.
        
    Returns:
        list[str]: Adapted metadata lines.
    """
    adapted_lines = []
    for line in meta_data:
        if line.startswith("DBREF"):
            adapted_line = handle_DBREF_line(line, region)
            adapted_lines.append(adapted_line)
        elif line.startswith("SEQRES"):
            parts = line.split()
            num_residues_region = region[1] - region[0] + 1
            parts[3] = str(num_residues_region)  # update number of residues
            adapted_line = " ".join(parts) + "\n"
            adapted_lines.append(adapted_line)
            # Note: Further adjustments to 3-letter code listing is not done yet.
        else:
            adapted_lines.append(line)
    return adapted_lines


def handle_DBREF_line(line: str, region: list[int]) -> str:
    """
    Adapts a DBREF line to reflect the extracted region.
    
    Args:
        line (str): Original DBREF line.
        region (list[int]): [start, end] residue indices.
        
    Returns:
        str: Adapted DBREF line.
    """
    region_start = region[0]
    region_end = region[1]
    parts = line.split()
    parts = line.split()
    pdb_start = int(parts[3])
    pdb_end = int(parts[4])
    uniprot_start = int(parts[8])
    uniprot_end = int(parts[9])
    offset_start = pdb_start - uniprot_start
    offset_end = pdb_end - uniprot_end
    uniprot_start_corrected = region_start - offset_start
    uniprot_end_corrected = region_end - offset_end
    parts[3] = str(region_start) # new pdb start reflects region start
    parts[4] = str(region_end) # new pdb end reflects region end
    parts[8] = str(uniprot_start_corrected) # new uniprot start reflects region start (corrected by possible offset)
    parts[9] = str(uniprot_end_corrected) # new uniprot end reflects region end (corrected by possible offset)
    adapted_line = " ".join(parts) + "\n"
    return adapted_line


def get_meta_data_of_structure_file(path_to_protein: str) -> list[str]:
    """
    Extracts the meta data from a PDB or CIF structure file. This includes header lines before the ATOM/HETATM lines.
    
    Args:
        path_to_protein (str): Path to the protein file.
        file_type (str): 'pdb' or 'cif'.
    Returns:
        list[str]: List of meta data lines.
    """
    with open(path_to_protein, "r") as f:
        lines = f.readlines()
        header_lines = []
        for line in lines:
            if line.startswith(("ATOM", "HETATM", "MODEL")):
                break
            header_lines.append(line)
    return header_lines


def download_fasta_record(uid: str) -> SeqRecord | bool:
    """
    Download a UniProt FASTA record by accession. If it is not available, try UniParc.

    Args:
        uid (str): UniProt accession.
    Returns:
        SeqRecord: downloaded FASTA record.
        bool: False if download failed.
    """
    # URLs for UniProt and UniParc
    uniprot_url = f"https://rest.uniprot.org/uniprotkb/{uid}.fasta"
    uniparc_url = f"https://rest.uniprot.org/uniparc/{uid}.fasta"
    for url in [uniprot_url, uniparc_url]:
        try:
            response = requests.get(url, timeout=10)
            if response.status_code == 200 and response.text.startswith(">"):
                fasta_io = StringIO(response.text)
                record = next(SeqIO.parse(fasta_io, "fasta"))
                return record
        except Exception as e:
            continue
    return False


def generate_fasta_from_uids_with_regions(uids_with_regions: dict, out_path: str, original_fasta=None):
    """
    Generates a FASTA file containing the sequences (cut down to their corresponding regions) of the given UIDs.
    The records will include the UID and region (if present) in the header.
    
    If original_fasta is provided, the method will extract sequences from it which are also part of the uids_with_regions.
    Otherwise the fasta file is generated from scratch and the Record IDs are of the format: UID or UID/start-end.

    Args:
        uids_with_regions (dict): Mapping {uid: [region_start, region_end] or None}.
        out_path (str): Path to the output FASTA file.
        original_fasta (str, optional): Path to an existing FASTA file with sequences.
    """
    # generate fasta with original fasta records
    if original_fasta is not None:
        uids = list(uids_with_regions.keys())
        copy_records_from_fasta(original_fasta, uids, out_path)
    # generate fasta with possible regions from scratch
    else:
        records = []
        for uid, region in uids_with_regions.items():
            downloaded_record = download_fasta_record(uid)
            if isinstance(downloaded_record, SeqRecord):
                record_with_region = add_region_to_record(downloaded_record, region)
                records.append(record_with_region)
            else:
                fallback_record = create_mok_up_record(uid, region)
                records.append(fallback_record)
        SeqIO.write(records, out_path, "fasta")


def create_mok_up_record(uid: str, region: list[int] | None) -> SeqRecord:
    """
    Creates a mock-up SeqRecord for a given UID and region.

    Args:
        uid (str): The UniProt ID.
        region (list[int] | None): The region to include in the ID.

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


def add_region_to_record(record: SeqRecord, region: list[int] | None) -> SeqRecord:
    """
    Takes a SeqRecord and modifies its ID to include the region if provided.

    Args:
        record (SeqRecord): The record to modify.
        region (list[int] | None): The region to include in the ID.

    Returns:
        SeqRecord: The modified record.
    """
    if region is None:
        region_str = ""
    else:
        region_str = f"/{int(region[0])}-{int(region[1])}"
    record.id = f"{record.id}{region_str}"
    record.description = ""
    return record


def copy_records_from_fasta(path_to_fasta: str, uids: list, out_path: str) -> str:
    """
    Copies the records with the given uids from the fasta file to a new fasta file.
    Args:
        path_to_fasta (str): Path to the input FASTA file.
        uids (list): List with UIDs to copy.
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
    