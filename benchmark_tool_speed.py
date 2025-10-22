import time
import os
from utils_for_PDB import fetch_pdbs
from StructSimComputer import StructSimComputer
from ToolType import ToolType
from InputFileType import InputFileType

"""
In this file the performance of different alignment tools of protein structures is benchmarked:
    - US-align
    - TM-align
    - Foldseek
This is done by running the tools on a predefined set of PDB files
and measuring the time taken for each tool with a timer.
"""

class Benchmark:
    def __init__(self, fasta_file=None, run_with_PDBs_for_benchmark=True):
        """
        Sets up the benchmark environment by downloading PDB files from a given FASTA file.
        If 'run_with_PDBs_for_benchmark' is False, the 'fasta_file' is mandatory, and it will overwrite the PDB files
        in the PDBs_for_benchmark directory. If 'run_with_PDBs_for_benchmark' is True, it will run the benchmark with already downloaded PDB files.
        In this case it does not matter if 'fasta_file' is given.
        """
        start = time.time()
        self.struct_sim_computer = StructSimComputer()
        self.results = {}
        if run_with_PDBs_for_benchmark:
            # run test with already downloaded PDB files
            if not os.path.exists("PDBs_for_benchmark"):
                print("PDBs_for_benchmark directory does not exist. Please rerun with a fasta file to download PDBs.")
                return
            else:
                self.data = "PDBs_for_benchmark"
                end = time.time()
                print(f"Benchmark environment setup complete after {end - start:.4f} seconds.")
            return
        else:
            if fasta_file is None:
                print("Please provide a fasta file for initial test.")
            else:
                print(f"Setting up benchmark environment with {fasta_file}...")
                fetch_pdbs(fasta_file, InputFileType.FASTA, "PDBs_for_benchmark")
                self.data = "PDBs_for_benchmark"
                end = time.time()
                print(f"Benchmark environment setup complete after {end - start:.4f} seconds. PDB files downloaded to './PDBs_for_benchmark'.")
    
    
    def run_benchmark(self):
        """
        Runs the benchmark for each tool in self.tools on all PDB files in self.data directory.
        """
        for tool in ToolType:
            print(f"Running benchmark for {tool.value}...")
            start = time.time()
            scores = self.struct_sim_computer.run(tool_type=tool, pdb_dir=self.data)
            end = time.time()
            total_time = end - start
            self.results[tool.name] = {
                'score': scores,
                'total_time': total_time
            }
            print(f"{tool.value} completed in {total_time:.4f} seconds")
            print(f"Scores: {scores}\n")
        print("Benchmark completed.")
        return self.results


# test
fasta_file = "./example_files/big_fasta_files/500.fasta"
benchmark = Benchmark(fasta_file=fasta_file, run_with_PDBs_for_benchmark=False)
results = benchmark.run_benchmark()
