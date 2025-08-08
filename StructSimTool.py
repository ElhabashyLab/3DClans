import subprocess
import os
import pandas as pd

"""
This module defines a base class for tools that can be used in a benchmarking context.
"""

class StructSimTool():
    def __init__(self, name: str, description: str, working_dir: str):
        self.name = name
        self.description = description
        self.working_dir = working_dir
        self.command = []
        self.output = None
        
        
    def start_run(self, pdb_dir):
        """
        Initializes the self.command list with the necessary parameters to run the tool and then returns _execute_run.
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
        Parses the output of the tool to extract the similarity scores.
        This method should be overridden by subclasses to implement specific parsing logic.
        """
        raise NotImplementedError("Subclasses should implement this method to parse output.")
        