from StructSimTool import StructSimTool
import os
import subprocess
"""
This class extends the Tool class to implement the USalign tool for protein structure comparison.
"""


class USalign(StructSimTool):
    def __init__(self):
        name = "USalign"
        description = "A tool for protein structure comparison using USalign."
        command = "USalign"
        self.flag_dir = "-dir" # specifies the directory containing PDB files
        self.flag_outfmt = "-outfmt 2" # specifies the output format
        super().__init__(name, description, command)
        
        
    def run(self, pdb_dir):
        """
        Run USalign on all PDB files in the directory pdb_dir.
        """
        # prepare files for USalign
        pdb_files = os.listdir(pdb_dir)
        with open("pdb_names.txt", 'w') as f:
            for pdb_name in pdb_files:
                f.write(f"{pdb_name}\n")
        # example command: USalign -dir PDBs_for_benchmark pdb_names.txt
        full_command = [self.command, self.flag_dir, pdb_dir, "pdb_names.txt", self.flag_outfmt]
        try:
            result = subprocess.run(full_command,
                                    capture_output=True,
                                    text=True,
                                    check=True)
            output = result.stdout
            # Parse similarity score from output
            score = self.parse_output(output)
        except subprocess.CalledProcessError as e:
            print(f"Error running {self.name} with {full_command}: {e}")
            score = None
        return score
    
    
    def parse_output(self, output):
        """
        Parses the output of the command to extract the similarity score.
        """
        return 1  # Placeholder for actual parsing logic, which should be implemented based on USalign's output format
    