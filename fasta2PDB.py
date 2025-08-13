from Bio import SeqIO
import requests
import os


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


def extract_uids_from_fasta(fasta_file):
    """
    Extract all UniProt IDs from Fasta file.
    All entries must contain a Uniprot_ID in the header.
    :return: a list which contains all Uniprot_IDs
    """
    uids = []
    for entry in SeqIO.parse(fasta_file, "fasta"):
        header = entry.id
        if "|" in header:
            uid = header.split("|")[1]
        else:
            # if no '|' is found in header, the header is assumed to be the uid
            uid = header
        uids.append(uid)
    return uids


def fetch_pdbs_from_uids(uniprot_ids, output_dir):
    """
    Fetches and stores PDB files in a specified output directory with a given list of uniprot_ids.
    """
    number_of_uids = len(uniprot_ids)
    number_of_failed_downloads = 0
    for uid in uniprot_ids:
        url = f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.pdb"
        if not download_file(url, output_dir + f"/{uid}.pdb"):
            number_of_failed_downloads += 1
    print(f"Downloaded {number_of_uids - number_of_failed_downloads} from {number_of_uids} PDB files successfully.")


"""
example test:
uids = extract_uids_from_fasta("example_files/combined.fasta")
fetch_pdbs_from_uids(uids, "./PDBs")
"""



