import subprocess
import os
import pandas as pd
from ClansFileGenerator import ClansFileGenerator
from InputFileType import InputFileType
from utils_for_structures_and_fasta import extract_uid_from_recordID
from ConfigFile import ConfigFile
        

def run_clans_headless(config_file: ConfigFile, path_clans_executable: str):
    """
    Runs clans executable in headless mode with a given config file.
    Args:
        config_file (ConfigFile): A config file object
        path_clans_executable (str): A path to the clans executable (jar file)
    Returns: str: A path to the output clans file as specified in the config file
    """
    arguments = config_file.read_config()
    outfile = arguments.get("saveto")
    arguments_of_conf = []
    for key, value in arguments.items():
        arguments_of_conf.append(f"-{key}")
        arguments_of_conf.append(str(value))
    print("\nrunning with " + " ".join(arguments_of_conf))
    args = ["java", "-Xmx4G", "-jar", path_clans_executable] + arguments_of_conf
    subprocess.run(args, check=True)
    return outfile


def generate_clans_file_seq_based(fasta_file_path: str, out_dir_path: str, blast_dir_path: str) -> str:
    """
    Generates a clans file based on sequence similarity from the given fasta file and blast results.
    _seq will be added to the output clans file name.
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
    clans_file_generator = ClansFileGenerator()
    out_path = os.path.join(out_dir_path, f"{input_file_name}_seq.clans")
    clans_file_path = clans_file_generator.generate_clans_file(scores_df, fasta_file_path, out_path)
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
                
                
def temp_add_pval_to_clans_file(clans_file_path: str, pval: float):
    """
    Adds pval to the parameter block in a CLANS file.
    If no <param> block exists, one is created.
    """
    # Lines you want to add (EDIT AS NEEDED)
    line_to_add = f"pval={pval}"

    with open(clans_file_path, "r") as f:
        lines = f.readlines()

    has_param_block = any(line.strip() == "<param>" for line in lines)

    if has_param_block:
        # Insert before </param>
        for i, line in enumerate(lines):
            if line.strip() == "</param>":
                lines.insert(i, line_to_add + "\n")
                break
    else:
        # Create param block after the first line (usually sequences=...)
        new_block = ["<param>\n"]
        new_block.append(line_to_add + "\n")
        new_block.append("</param>\n")

        # Insert after header line
        lines[1:1] = new_block

    with open(clans_file_path, "w") as f:
        f.writelines(lines)
    