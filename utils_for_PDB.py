import re
from Bio import SeqIO, Entrez
import Bio.PDB as bpdb
Entrez.email = "aron.wichtner@tuebingen.mpg.de"
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
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        response = requests.get(url, stream=True, timeout=10)
        # raise exception if request failed
        response.raise_for_status()
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
    It returns a dict containing uid : [region_start, region_end] of the sequences that have been successfully downloaded.
    
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
    output_path = os.path.join(output_dir, f"{uid}.{file_format}")
    success = download_file(url, output_path)
    if success and region:
        ...
    return success


class ResidueSelect(Select):
    def __init__(self, start, end):
        self.start = start
        self.end = end
    def accept_residue(self, residue):
        res_id = residue.get_id()[1]  # residue number
        return self.start <= res_id <= self.end


def extract_region_of_protein(path_to_protein: str, file_type: str, region: list[int]):
    """
    Cuts the protein to a desired start and end Aminoacid as specifies in region.
    The protein can be of file type PDB or cif.

    Args:
        path_to_protein (str): Path to the file of the protein to be cut.
        file_type (str): 'pdb' or 'cif'
        region (list[int]): Residue range [start, end]
    """
    file_name 
    parser = PDBParser()
    structure = parser.get_structure("AF", f"{uniprot_id}.pdb")

    io = PDBIO()
    io.set_structure(structure)
    io.save(f"{uniprot_id}_region.pdb", ResidueSelect(10, 30))  # residues 10-30
    


def download_fasta_record(uid: str):
    """
    Download protein fasta record with given uid
    
    Args: 
        uid (str): Protein uid.
    Returns:
        SeqRecord: downloaded record.
    """
    with Entrez.efetch(db="protein", id=uid, rettype="fasta", retmode="text") as handle:
        record = SeqIO.read(handle, "fasta")
    return record
        

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
    records = []

    # generate fasta with original fasta records
    if original_fasta is not None:
        original_records = list(SeqIO.parse(original_fasta, "fasta"))
        for uid, region in uids_with_regions.items():
            for record in original_records:
                if uid in record.id:
                    records.append(record)
    # generate fasta from scratch
    else:
        for uid, region in uids_with_regions.items():
            record = download_fasta_record(uid)
            records.append(record)
    # write records to fasta file
    SeqIO.write(records, out_path, "fasta")
    