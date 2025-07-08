from Tool import Tool
import os
import subprocess
"""
This class extends the Tool class to implement the TMalign tool for protein structure alignment.
"""


class TMalign(Tool):
    def __init__(self, name, description, command):
        super().__init__(name, description, command)

    def run(self, pdb_dir="PDBs_for_benchmark"):
        """
        Runs the TMalign tool on all PDB files in the specified directory.
        """
        filesnames = os.listdir(pdb_dir)
        for pdbi in filesnames:
            for pdbj in filesnames:
                if pdbi != pdbj:
                    pdb_file_path_i = os.path.join(pdb_dir, pdbi)
                    pdb_file_path_j = os.path.join(pdb_dir, pdbj)
                    full_command = f"{self.command} {pdb_file_path_i} {pdb_file_path_j}"
                    subprocess.run(full_command, shell=True)
