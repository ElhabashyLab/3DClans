import cProfile
import os
from fasta2PDB import *
from TMalign import TMalign
"""
In this file/directory the performance of different alignment tools of protein structures is benchmarked:
    - US-align
    - TM-align
    - Foldseek
    - 3D Blast
    - Graph Align
This is done by running the tools on a predefined set of PDB files
and measuring the time taken for each tool with cProfile.
"""


def setup_benchmark(fasta_file=None, is_initial_test=True):
    """
    Sets up the benchmark environment by downloading PDB files from a given FASTA file.
    If 'is_initial_test' is True, the 'fasta_file' is mandatory, and it will download the PDB files
    and delete the content of the PDB directory if it already exists.
    """
    tools = set_up_tools()
    if not is_initial_test:
        # run test with already downloaded PDB files
        return
    else:
        if fasta_file is None:
            print("please provide a fasta file for initial test")
        else:
            print(f"Setting up benchmark environment with {fasta_file}...")
            if os.path.exists("PDBs_for_benchmark"):
                delete_dir_content("PDBs_for_benchmark")
            else:
                os.makedirs("PDBs_for_benchmark")
            UIDS = extract_uids_from_fasta(fasta_file)
            fetch_pdbs_from_uids(UIDS, "PDBs_for_benchmark")
            print("Benchmark environment setup complete. PDB files downloaded to './PDBs_for_benchmark'.")


def set_up_tools():
    """
    Sets up the tools required for benchmarking.
    """
    tools = []
    # set up each tool by extending the Tool class
    tmalign = TMalign("TM-align",
                      "Tool for protein structure alignment using TM-score",
                      "$TMalign")
    tools.append(tmalign)
    return tools


def run_benchmark_tool(tool_name, pdb_dir="PDBs_for_benchmark"):
    """
    Runs the specified benchmark tool on all PDB files in the given directory and profiles its performance.
    :param tool_name: Name of the tool to benchmark (e.g., 'US-align', 'TM-align', etc.)
    :param pdb_dir: Directory containing PDB files to benchmark
    """
    print(f"Running {tool_name} benchmark...")
    cProfile.run(f"{tool_name}.run('{pdb_dir}')", f"{tool_name}_benchmark.prof")
    print(f"{tool_name} benchmark completed.")


def delete_dir_content(dir_path):
    """
    Deletes the content of the specified directory.
    """
    if os.path.exists(dir_path):
        for file in os.listdir(dir_path):
            file_path = os.path.join(dir_path, file)
            try:
                if os.path.isfile(file_path) or os.path.islink(file_path):
                    os.unlink(file_path)
                elif os.path.isdir(file_path):
                    os.rmdir(file_path)
            except Exception as e:
                print(f"Failed to delete {file_path}. Reason: {e}")
    else:
        print(f"Directory {dir_path} does not exist.")


INPUT_FILE = "./input/small_dataset.a2m"

setup_benchmark(INPUT_FILE)
