import subprocess
import os
"""
This module defines a base class for tools that can be used in a benchmarking context.
"""


class Tool:
    def __init__(self, name: str, description: str, command: str):
        self.name = name
        self.description = description
        self.command = command
        self.output = None

    def run_all_against_all(self, pdb_dir="PDBs_for_benchmark"):
        """
        Runs the Tool on all PDB files in the specified directory.
        """
        TM_scores = []
        filesnames = os.listdir(pdb_dir)
        for pdbi in filesnames:
            for pdbj in filesnames:
                if pdbi != pdbj:
                    pdb_file1_path = os.path.join(pdb_dir, pdbi)
                    pdb_file2_path = os.path.join(pdb_dir, pdbj)
                    tm_score = self.run(pdb_file1_path, pdb_file2_path)
                    TM_scores.append(tm_score)
        return TM_scores

    def run(self, pdb_file1_path, pdb_file2_path):
        """
        Runs the Tool on two specified PDB files and return TM-score.
        """
        full_command = f"{self.command} {pdb_file1_path} {pdb_file2_path}"
        try:
            result = subprocess.run(full_command,
                                    capture_output=True,
                                    text=True,
                                    check=True)
            output = result.stdout
            # Parse TM-score from output
            tm_score = None
        except subprocess.CalledProcessError as e:
            print(f"Error running {self.name} with {full_command}: {e}")
            tm_score = None

        return tm_score




