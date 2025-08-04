import subprocess
import os
import pandas as pd

"""
This module defines a base class for tools that can be used in a benchmarking context.
"""

class StructSimTool():
    def __init__(self, name: str, description: str):
        self.name = name
        self.description = description
        self.command = []
        self.output = None
        
        
    def start_run(self, pdb_dir):
        """
        RuInitializes the tool to run on the specified PDB directory.
        This method should be overridden by subclasses to implement specific tool logic.
        """
        raise NotImplementedError("Subclasses should implement this method to run the tool.")
    
    
    def _execute_run(self):
        """
        Executes self.command in a subprocess and captures the output.
        """
        try:
            result = subprocess.run(self.command,
                                    capture_output=True,
                                    text=True,
                                    check=True)
            self.output = result.stdout
            scores = self._parse_output()
        except subprocess.CalledProcessError as e:
            print(f"Error running {self.name} with {self.command}: {e}")
            scores = None
        return scores
    
    
    def _parse_output(self):
        """
        Parses the output of the command to extract the similarity score.
        This method should be overridden by subclasses to implement specific parsing logic.
        """
        raise NotImplementedError("Subclasses should implement this method to parse output.")
        
"""
    def run_all_against_all(self, pdb_dir):
        
        Runs the Tool on all PDB files in the specified directory to compute pairwise similarity scores.
        Returns a data frame with similarity scores for each pair of PDB files.
        
        scores = []
        filenames = os.listdir(pdb_dir)
        for i in range(len(filenames)):
            for j in range(i + 1, len(filenames)):
                if i != j:
                    pdb_file1_path = os.path.join(pdb_dir, filenames[i])
                    pdb_file2_path = os.path.join(pdb_dir, filenames[j])
                    pair_score = self.run(pdb_file1_path, pdb_file2_path)
                    scores.append({
                        'file1': filenames[i],
                        'file2': filenames[j],
                        'score': pair_score
                    })
        # Convert scores to a DataFrame
        scores_df = pd.DataFrame(scores)
        return scores
"""
    