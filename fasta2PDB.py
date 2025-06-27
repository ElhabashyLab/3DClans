from Bio import SeqIO
import requests


def extract_UIDS_from_fasta(fasta_file):
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
            uid = header
        uids.append(uid)
    return uids


def download_file(url, output_path):
    """
    Download a file from url to a output_path.
    """
    try:
        response = requests.get(url, stream=True, timeout=10)
        response.raise_for_status()
        with open(output_path, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                f.write(chunk)
        return True
    except Exception as e:
        print(f"Failed to download {url}: {str(e)}")
        return False


def fetch_pdbs(uniprot_ids, output_dir):
    urls = [f"https://alphafold.ebi.ac.uk/files/AF-{uid}-F1-model_v4.pdb" for uid in uniprot_ids]
    for url in urls:
        download_file(url, ...)



print(extract_UIDS_from_fasta("example_files/example.fasta"))