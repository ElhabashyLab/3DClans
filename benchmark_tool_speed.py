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
This is done by running the code of `benchmark_tool_speed.py` on a predefined set of protein structures in a folder
or by supplying a valid input-file. The Benchmark measures and records the time taken for each tool.
"""

class Benchmark:
    def __init__(self, input_file=None, input_file_type=None, run_with_structures_for_benchmark=True):
        """
        Sets up the benchmark environment by downloading protein structure files from a given input-file.
        If 'run_with_structures_for_benchmark' is False, the 'input_file' is mandatory, and it will overwrite the protein structure files
        in the structures_for_benchmark directory. In this case it is also mandatory to specify the 'input_file_type'. 
        If 'run_with_structures_for_benchmark' is True, it will run the benchmark with already downloaded protein structure files.
        In this case it does not matter if 'input_file' or 'input_file_type' is given.
        """
        start = time.time()
        self.struct_sim_computer = StructSimComputer()
        self.results = {}
        if run_with_structures_for_benchmark:
            # run test with already downloaded structure files
            if not os.path.exists("structures_for_benchmark"):
                print("structures_for_benchmark directory does not exist. Please rerun with a fasta file to download protein structure files.")
                return
            else:
                self.data = "structures_for_benchmark"
                end = time.time()
                print(f"Benchmark environment setup complete after {end - start:.4f} seconds.")
            return
        else:
            if input_file is None or input_file_type is None:
                print("Please provide a fasta file with a correct input_file_type for initial test.")
            else:
                print(f"Setting up benchmark environment with {input_file}...")
                uids_with_regions = fetch_pdbs(input_file, input_file_type, "structures_for_benchmark")
                self.data = "structures_for_benchmark"
                end = time.time()
                print(f"Benchmark environment setup complete after {end - start:.4f} seconds. Structure files downloaded to './structures_for_benchmark'.")
    
    
    def run_benchmark(self):
        """
        Runs the benchmark for each tool in self.tools on all structure files in self.data directory.
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
input_file = "./example_files/small_tsv_files/test_1.tsv"
benchmark = Benchmark(input_file=input_file, input_file_type=InputFileType.TSV, run_with_structures_for_benchmark=False)
results = benchmark.run_benchmark()
