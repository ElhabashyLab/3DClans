from Bio import SeqIO
import requests
import os
import shutil


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
        Deletes the content of the specified directory.
        """
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
    Cleans the fasta file by removing entries that do not have a corresponding PDB file.
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


def fetch_pdbs(fasta_file, output_dir):
    """
    Fetches and stores PDB files in a specified output directory with a given fasta file.
    The contents of the output dir will be overwritten.
    It also returns a cleaned fasta file containing only the sequences that have been successfully downloaded.
    """
    parsed_fasta = SeqIO.parse(fasta_file, "fasta")
    number_of_uids = 0
    number_of_failed_downloads = 0
    # get dir of fasta_file
    fasta_dir = os.path.dirname(fasta_file)
    cleaned_fasta_file = os.path.join(fasta_dir, "cleaned.fasta")
    # delete old PDB files if they exist
    if os.path.exists(output_dir):
        delete_dir_content(output_dir)
    else:
        os.makedirs(output_dir)
    with open(cleaned_fasta_file, 'w') as file:
        for record in parsed_fasta:
            number_of_uids += 1
            uid = extract_uid_from_recordID(record.id)
            url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.pdb"
            if download_file(url, output_dir + f"/{uid}.pdb"):
                file.write(f">{record.id}\n{str(record.seq)}\n")
            else:
                number_of_failed_downloads += 1
    print(f"Downloaded {number_of_uids - number_of_failed_downloads} from {number_of_uids} PDB files successfully.")
    return cleaned_fasta_file


"""
example test:
uids = extract_uids_from_fasta("example_files/combined.fasta")
fetch_pdbs_from_uids(uids, "./PDBs")
"""



