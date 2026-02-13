import os
import requests
import pandas as pd
from Bio import SeqIO
from Bio.PDB.PDBIO import PDBIO
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
from Bio.PDB.PDBIO import Select
from clans3d.core.input_file_type import InputFileType
from clans3d.utils.file_utils import download_file, reset_dir_content
from clans3d.utils.fasta_utils import extract_uid_from_recordID, extract_region_from_record


def fetch_pdbs(input_file_path: str, input_file_type: InputFileType, output_dir: str) -> dict[str, tuple[int, int] | None]:
    """
    Fetches and stores PDB files, specified in a given input_file in a output directory.
    The contents of the output dir will be overwritten.
    It returns a dict containing {uid : (region_start, region_end)} of the sequences that have been successfully downloaded.
    
    Args:
        input_file_path: Path to the input file.
        input_file_type: Specifies the type of the input file.
        output_dir: The directory where the downloaded PDBs are stored.
        
    Returns:
        dict (str, tuple[int, int] | None): Containing downloaded uids together with their regions.
    """    
    reset_dir_content(output_dir)
    print(f"Downloading structure files in \"{output_dir}\"...")
    if input_file_type is InputFileType.FASTA or input_file_type is InputFileType.A2M:
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
    Downloads the structures of the uids in the given tsv file.
    It also creates a dictionary with an entry for each downloaded record. {uid : [region]}
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
            region = (int(start), int(end))
        if download_alphafold_structure(uid, output_dir, region):
            successful_downloads += 1
            uids_with_regions[uid] = region
    print(f"Downloaded {successful_downloads} from {total_uids} PDB files successfully.")
    return uids_with_regions


def download_alphafold_structure(
    id: str,
    output_dir: str,
    region: tuple[int, int] | None,
    file_format: str = "pdb",
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
        file_format (str): Structure format to download ("pdb" or "cif").

    Returns:
        bool: True if the structure was successfully downloaded (and processed),
              False otherwise.
    """
    if file_format not in {"pdb", "cif"}:
        raise ValueError("file_format must be 'pdb' or 'cif'")
    
    api_url = f"https://alphafold.ebi.ac.uk/api/prediction/{id}"
    try:
        response = requests.get(api_url, timeout=30)
        response.raise_for_status()
        predictions = response.json()
    except Exception:
        print(f"Failed to download {id} from {api_url}")
        return False
    if not predictions:
        print(f"Failed to download {id} from {api_url}")
        return False
    
    model_info = predictions[0]
    if file_format == "pdb":
        structure_url = model_info.get("pdbUrl")
    else:
        structure_url = model_info.get("cifUrl")
    if not structure_url:
        print(f"Failed to download {id} from {api_url}")
        return False

    path_to_structure = os.path.join(output_dir, f"{id}.{file_format}")
    success = download_file(structure_url, path_to_structure)
    if not success:
        print(f"Failed to download {id} from {api_url}")
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
    without corrupting metadata. Writes output to a new file by default.
    
    Args:
        path_to_protein (str): Path to the protein file.
        file_type (str): 'pdb' or 'cif'.
        region (tuple[int, int]): (start, end) residue indices to extract.
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
        out.writelines(meta_data_corrected)
        io = PDBIO()
        io.set_structure(structure)
        io.save(out, select=SelectRegion())
    return output_path


def adapt_metadata_to_region(meta_data: list[str], region: tuple[int, int]) -> list[str]:
    """
    Adapts the metadata lines of a structure file to reflect the extracted region.
    
    Args:
        meta_data (list[str]): Original metadata lines.
        region (tuple[int, int]): [start, end] residue indices.
        
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


def handle_DBREF_line(line: str, region: tuple[int, int]) -> str:
    """
    Adapts a DBREF line to reflect the extracted region.
    
    Args:
        line (str): Original DBREF line.
        region (tuple[int, int]): [start, end] residue indices.
        
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
