from Bio import SeqIO
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


def fetch_pdbs(input_file_path: str, input_file_type: InputFileType, output_dir: str) -> str:
    """
    Fetches and stores PDB files, specified in a given input_file in a output directory.
    The contents of the output dir will be overwritten.
    It also returns a cleaned input_file containing only the sequences that have been successfully downloaded and saves it in the same dir as the given input_file.
    
    Args:
        input_file_path: Path to the input file.
        input_file_type: Specifies the type of the input file.
        output_dir: The directory where the downloaded PDBs are stored.
        
    Returns:
        cleaned_input_file_path: Path to cleaned version of the input_file containing only successfully downloaded elements.
    """    
    input_file_name = os.path.basename(input_file_path).split(".")[0]
    input_file_dir = os.path.dirname(input_file_path)
    cleaned_input_file_path = os.path.join(input_file_dir, f"{input_file_name}_cleaned.{input_file_type.value}")
    delete_dir_content(output_dir)    
    # Process based on file type
    if input_file_type is InputFileType.FASTA:
        return _process_fasta_file(input_file_path, cleaned_input_file_path, output_dir)
    elif input_file_type is InputFileType.TSV:
        return _process_tsv_file(input_file_path, cleaned_input_file_path, output_dir)
    else:
        raise ValueError(f"Unsupported input file type: {input_file_type}")


def _process_fasta_file(input_file_path: str, cleaned_file_path: str, output_dir: str) -> str:
    """Process FASTA file and download corresponding PDB files."""
    successful_downloads = 0
    total_records = 0
    with open(cleaned_file_path, 'w') as cleaned_file:
        for record in SeqIO.parse(input_file_path, "fasta"):
            total_records += 1
            uid = extract_uid_from_recordID(record.id)
            if _download_alphafold_structure(uid, output_dir):
                SeqIO.write(record, cleaned_file, "fasta")
                successful_downloads += 1
    print(f"Downloaded {successful_downloads} from {total_records} PDB files successfully.")
    return cleaned_file_path


def _process_tsv_file(input_file_path: str, cleaned_file_path: str, output_dir: str) -> str:
    """Process tsv file and download corresponding PDB files."""
    df = pd.read_csv(input_file_path, sep='\t')
    successful_downloads = 0
    total_uids = len(df)
    def download_and_check(row):
        nonlocal successful_downloads
        uid = row['uid']
        if _download_alphafold_structure(uid, output_dir):
            successful_downloads += 1
            return True
        return False
    successful_mask = df.apply(download_and_check, axis=1)
    successful_df = df[successful_mask]
    successful_df.to_csv(cleaned_file_path, sep='\t', index=False)
    print(f"Downloaded {successful_downloads} from {total_uids} PDB files successfully.")
    return cleaned_file_path


def _download_alphafold_structure(uid: str, output_dir: str) -> bool:
    """Download AlphaFold structure for a given UID."""
    url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v6.pdb" # This url might change over time depending on model
    output_path = os.path.join(output_dir, f"{uid}.pdb")
    return download_file(url, output_path)
