import cProfile
import os
from fasta2PDB import extract_uids_from_fasta, fetch_pdbs_from_uids
from USalign import USalign
from TMalign import TMalign
from Foldseek import Foldseek
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

class Benchmark:
    def __init__(self, fasta_file=None, run_with_PDBs_for_benchmark=True):
        """
        Sets up the benchmark environment by downloading PDB files from a given FASTA file.
        If 'run_with_PDBs_for_benchmark' is False, the 'fasta_file' is mandatory, and it will overwrite the PDB files
        in the PDBs_for_benchmark directory. If 'run_with_PDBs_for_benchmark' is True, it will run the benchmark with already downloaded PDB files.
        In this case it does not matter if 'fasta_file' is given.
        """
        self.tools = self._set_up_tools()
        self.results = {}
        if run_with_PDBs_for_benchmark:
            # run test with already downloaded PDB files
            if not os.path.exists("PDBs_for_benchmark"):
                print("PDBs_for_benchmark directory does not exist. Please rerun with a fasta file to download PDBs.")
                return
            else:
                self.data = "PDBs_for_benchmark"
            return
        else:
            if fasta_file is None:
                print("Please provide a fasta file for initial test.")
            else:
                print(f"Setting up benchmark environment with {fasta_file}...")
                if os.path.exists("PDBs_for_benchmark"):
                    self._delete_dir_content("PDBs_for_benchmark")
                else:
                    os.makedirs("PDBs_for_benchmark")
                UIDS = extract_uids_from_fasta(fasta_file)
                fetch_pdbs_from_uids(UIDS, "PDBs_for_benchmark")
                self.data = "PDBs_for_benchmark"
                print("Benchmark environment setup complete. PDB files downloaded to './PDBs_for_benchmark'.")
    
    
    def run_benchmark(self):
        """
        Runs the benchmark for each tool in self.tools on all PDB files in self.data directory.
        """
        for tool in self.tools:
            print(f"Running benchmark for {tool.name}...")
            profiler = cProfile.Profile()
            profiler.enable()
            scores = tool.start_run(self.data)
            profiler.disable()
            total_time = sum(stat.totaltime for stat in profiler.getstats())
            self.results[tool.name] = {
                'score': scores,
                'total_time': total_time
            }
            print(f"{tool.name} completed in {total_time:.4f} seconds")
        return self.results


    def _set_up_tools(self):
        """
        Sets up the tools required for benchmarking.
        """
        tools = [foldseek := Foldseek(), usalign := USalign(), tmalign := TMalign()]
        return tools


    def _delete_dir_content(self, dir_path):
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


# test
fasta_file = "./example_files/small_fasta_files/small_dataset.a2m"
benchmark = Benchmark(fasta_file=fasta_file, run_with_PDBs_for_benchmark=False)
results = benchmark.run_benchmark()
