import subprocess
import os
import pandas as pd
from ClansFileGenerator import ClansFileGenerator
from InputFileType import InputFileType
from utils_for_structures_and_fasta import extract_uid_from_recordID


def run_clans_headless_from_c_file(input_file: str, output_file: str, recovered_clans_path: str, rounds: int = 100):
    """
    Runs recovered clans in headless mode on the given clans file and saves the output to the output file.
    Args:
        input_file: A path to the input clans file
        output_file: A path to the output clans file
        recovered_clans_path: A path to the recovered clans project directory
        rounds: The number of rounds to run clans for.
    Returns: None
    """
    gradlew = "./gradlew"
    args = [
        gradlew, "run", "--no-daemon",
        f'--args=-nographics T -load {input_file} -dorounds {rounds} -saveto {output_file}'
    ]
    subprocess.run(args, cwd=recovered_clans_path, check=True)


def run_clans_headless_from_f_file(input_file: str, output_file: str, recovered_clans_path: str, blast_dir: str, rounds: int = 100):
    """
    Creates a clans file from the given fasta file using sequence similarity,
    and then runs recovered clans in headless mode on that clans file and saves the output to the output file.
    Default settings for blast will be used.
    Args:
        input_file (str): A path to the input fasta file
        output_file (str): A path to the output clans file
        recovered_clans_path (str): A path to the recovered clans project directory
        blast_dir (str): A path to the directory where the blast database and results will be stored
        rounds (int): The number of rounds to run clans for.
    Returns: None
    """
    out_dir_path = os.path.dirname(output_file)
    seq_based_clans_file = generate_clans_file_seq_based(input_file, out_dir_path, blast_dir) 
    run_clans_headless_from_c_file(seq_based_clans_file, output_file, recovered_clans_path, rounds)


def generate_clans_file_seq_based(fasta_file_path: str, out_dir_path: str, blast_dir_path: str) -> str:
    """
    Generates a clans file based on sequence similarity from the given fasta file and blast results.
    The clans files is generated using the default parameters of the CLANS-web-tool.
    Args:
        fasta_file_path (str): A path to the input fasta file
        out_dir_path (str): A path to the output directory where the clans file will be saved
        blast_dir_path (str): A path to the directory where the blast database and results will be stored
    Returns:
        str: A path to the generated clans file
    """
    input_file_name = os.path.splitext(os.path.basename(fasta_file_path))[0]
    blast_results_path = os.path.join(blast_dir_path, f"{input_file_name}_blast.tsv")
    [blast_results_path, outformat] = blast_fasta(fasta_file_path, blast_results_path, blast_dir_path)
    blast_results_df = pd.read_csv(blast_results_path, sep="\t", names=["PDBchain1", "PDBchain2", "score"])
    blast_results_df["PDBchain1"] = blast_results_df["PDBchain1"].apply(extract_uid_from_recordID)
    blast_results_df["PDBchain2"] = blast_results_df["PDBchain2"].apply(extract_uid_from_recordID)
    blast_results_df = blast_results_df[blast_results_df['PDBchain1'] != blast_results_df['PDBchain2']] # remove self-hits
    scores_df = blast_results_df.drop_duplicates(subset=["PDBchain1", "PDBchain2"]).reset_index(drop=True)
    clans_file_generator = ClansFileGenerator(out_dir_path)
    clans_file_path = clans_file_generator.generate_clans_file(scores_df, fasta_file_path)
    return clans_file_path


def blast_fasta(fasta_file: str, output_file: str, working_dir: str):
    """
    Runs blast on the given fasta file and saves the output to the output file.
    Args:
        fasta_file (str): A path to the input fasta file
        output_file (str): A path to the output blast file
        working_dir (str): A path to the working directory where the blast database will be created
    Returns: None
    """
    tmp_db = os.path.join(working_dir, "temp_blast_db") 
    outfmt = "6 qseqid sseqid evalue"
    subprocess.run([ "makeblastdb", "-in", fasta_file, "-dbtype", "prot", "-out", tmp_db ], check=True)
    blast_cmd = [ "blastp",
                 "-query", fasta_file,
                 "-db", tmp_db,
                 "-out", output_file,
                 "-outfmt", outfmt,
                 "-evalue", "1e-4",
                 "-matrix", "BLOSUM62"] 
    subprocess.run(blast_cmd, check=True)
    return [output_file, outfmt.split(" ")[1:]]
    

def run_clans_headless(input_output_files: dict, input_file_type: InputFileType, recovered_clans_path: str, rounds: int = 100000, blast_dir = None):
    """
    Runs recovered clans in headless mode on each of the given input files and saves the output to the corresponding output files.
    The input_output_files dict should contain input_file_path: output_file_path pairs.
    Args:
        input_output_files (dict): A dict containing input_file_path: output_file_path pairs
        input_file_type (InputFileType): The type of the input files. Can be either InputFileType.CLANS or InputFileType.FASTA.
        recovered_clans_path (str): A path to the recovered clans project directory
        rounds (int): The number of rounds to run clans for.
        blast_dir (str): A path to the directory where the blast database and results will be stored. Only needed if input_file_type is "fasta".
        clans_generator (ClansFileGenerator): An optional ClansFileGenerator object to use for generating clans files from fasta files. Must be provided if input_file_type is "fasta".
    Returns: None
    """
    for input_file, output_file in input_output_files.items():
        if input_file_type == InputFileType.CLANS:
            run_clans_headless_from_c_file(input_file, output_file, recovered_clans_path, rounds)
        elif input_file_type == InputFileType.FASTA:
            if blast_dir is None:
                raise ValueError("blast_dir must be provided when input_file_type is 'fasta'")
            run_clans_headless_from_f_file(input_file, output_file, recovered_clans_path, blast_dir, rounds)
        else:
            raise ValueError(f"Invalid input file type: {input_file_type}. Supported types are {InputFileType.CLANS} and {InputFileType.FASTA}.")
        _remove_colorcutoffs_colorarr(output_file)


def _remove_colorcutoffs_colorarr(clans_file):
        """
        Fixes a temporary bug in the recovered clans code.
        This function is a temporary workaround and should be removed once the bug is fixed in the recovered clans code.
        After the clustering process the clans files contain 2 lines starting with colorcutoffs and colorarr with missing values.
        This function removes those lines from the clans files in the given directory so they can be loaded again.
        Args:
            clans_files: A path to the the clans file to be fixed.
        """
        with open(clans_file, "r") as f:
            lines = f.readlines()
        with open(clans_file, "w") as f:
            for line in lines:
                if line.startswith("colorcutoffs") or line.startswith("colorarr"):
                    continue # remove buggy lines
                f.write(line)
    

#test
#PATH_TO_RECOVERED_CLANS = "/home/aronw/Development/clans-recovered"
#input_file = "/home/aronw/Development/clans_files_seqsim/dataset_1_cleaned.clans"
#output_file = "/home/aronw/Development/clans_files_seqsim/dataset_1_cleaned_updated.clans"
#dict_test = {input_file: output_file}
#run_clans_headless(PATH_TO_RECOVERED_CLANS, dict_test)
