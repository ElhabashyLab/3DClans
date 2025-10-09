import subprocess
import os
import pandas as pd
from ClansFileGenerator import ClansFileGenerator


# This is the path to the recovered clans project directory (github-repo: https://github.com/AronWichtner/clans-recovered.git)


def run_clans_headless_from_c_file(recovered_clans_path: str, input_file: str, output_file: str, rounds: int = 100):
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


def run_clans_headless_from_f_file(recovered_clans_path: str, input_file: str, output_file: str, blast_dir: str, clans_file_generator: ClansFileGenerator, rounds: int = 100):
    """
    Runs recovered clans in headless mode on the given fasta file and saves the output to the output file.
    From the fasta file, a clans file will be created using blast, and then clans will be run on that clans file.
    Default settings for blast will be used.
    Args:
        input_file: A path to the input fasta file
        output_file: A path to the output clans file
        recovered_clans_path: A path to the recovered clans project directory
        blast_dir: A path to the directory where the blast database and results will be stored
        rounds: The number of rounds to run clans for.
    Returns: None
    """
    # get pairwise blast results from fasta file
    input_file_name = os.path.splitext(os.path.basename(input_file))[0]
    blast_results_path = os.path.join(blast_dir, f"{input_file_name}_blast.tsv")
    [blast_results_path, outformat] = blast_fasta(input_file, blast_results_path, blast_dir)
    blast_results_df = pd.read_csv(blast_results_path, sep="\t", names=["PDBchain1", "PDBchain2", "score"])
    blast_results_df = blast_results_df[blast_results_df['PDBchain1'] != blast_results_df['PDBchain2']] # remove self-hits
    blast_results_df = blast_results_df.drop_duplicates(subset=["PDBchain1", "PDBchain2"]).reset_index(drop=True)
    clans_file_path = clans_file_generator.generate_clans_file(blast_results_df, input_file)
    run_clans_headless_from_c_file(recovered_clans_path, clans_file_path, output_file, rounds)


def blast_fasta(fasta_file: str, output_file: str, working_dir: str):
    """
    Runs blast on the given fasta file and saves the output to the output file.
    Args:
        fasta_file: A path to the input fasta file
        output_file: A path to the output blast file
        working_dir: A path to the working directory where the blast database will be created
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
    

def run_clans_headless(recovered_clans_path: str, input_output_files: dict, input_file_type: str = "clans", rounds: int = 100000, blast_dir = None, clans_generator = None):
    """
    Runs recovered clans in headless mode on each of the given input files and saves the output to the corresponding output files.
    The input_output_files dict should contain input_file: output_file pairs.
    Args:
        input_output_files: A dict containing input_file: output_file pairs
        recovered_clans_path: A path to the recovered clans project directory
        input_file_type: The type of the input files. Can be either "clans" or "fasta".
        rounds: The number of rounds to run clans for.
        blast_dir: A path to the directory where the blast database and results will be stored. Only needed if input_file_type is "fasta".
        clans_generator: An optional ClansFileGenerator object to use for generating clans files from fasta files. Must be provided if input_file_type is "fasta".
    Returns: None
    """
    for input_file, output_file in input_output_files.items():
        if input_file_type == "clans":
            run_clans_headless_from_c_file(recovered_clans_path, input_file, output_file, rounds)
        elif input_file_type == "fasta":
            if blast_dir is None:
                raise ValueError("blast_dir must be provided when input_file_type is 'fasta'")
            if clans_generator is None:
                raise ValueError("clans_generator must be provided when input_file_type is 'fasta'")
            run_clans_headless_from_f_file(recovered_clans_path, input_file, output_file, blast_dir, clans_generator, rounds)
        else:
            raise ValueError(f"Invalid input file type: {input_file_type}. Supported types are 'clans' and 'fasta'.")


#test
PATH_TO_RECOVERED_CLANS = "/home/aronw/Development/clans-recovered"
input_file = "/home/aronw/Development/Clans-3D/Scores_Evaluation/clans_files_seqsim/dataset_1_cleaned.clans"
output_file = "/home/aronw/Development/Clans-3D/Scores_Evaluation/clans_files_seqsim/dataset_1_cleaned_updated.clans"
dict_test = {input_file: output_file}
#run_clans_headless(PATH_TO_RECOVERED_CLANS, dict_test)
