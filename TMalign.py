from StructSimTool import StructSimTool
import os
import subprocess
import pandas as pd
from io import StringIO
"""
This class extends the StructSimTool class to implement the TMalign tool for protein structure comparison.
"""


class TMalign(StructSimTool):
    def __init__(self):
        description = "A tool for protein structure comparison using TMalign."
        super().__init__("TMalign", description)
        self.flag_dir = "-dir" # specifies the directory containing PDB files
        self.flag_outfmt = "-outfmt" # specifies the output format
        self.outfmt = "2"
        
        
    def run(self, pdb_dir):
        """
        Runs the tool on the specified PDB directory.
        """
        # prepare files for TMalign
        pdb_files = os.listdir(pdb_dir)
        with open("pdb_names.txt", 'w') as f:
            for pdb_name in pdb_files:
                f.write(f"{pdb_name}\n")
        # example command: "TMalign -dir PDBs_for_benchmark/ pdb_names.txt -outfmt 2"
        self.command = [self.name, self.flag_dir, pdb_dir + "/", "pdb_names.txt", self.flag_outfmt, self.outfmt]
        return self._execute_run()

    
    def _parse_output(self):
        """
        Parses the output of the command to extract the similarity score.
        """
        try:
            df = pd.read_csv(StringIO(self.output), sep='\t')
            if not df.empty:
                df.columns = [col.lstrip('#') for col in df.columns]
                scores = df[['PDBchain1', 'PDBchain2', 'TM1', 'TM2']]
                return scores
            else:
                print("Failed to parse output: DataFrame is empty.")
                return None
        except Exception as e:
            print(f"Error parsing output: {e}")
            return False
        