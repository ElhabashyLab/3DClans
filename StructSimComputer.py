from USalign import USalign
from TMalign import TMalign
from Foldseek import Foldseek
from ToolType import ToolType
import os


class StructSimComputer:
    """
    This class is responsible for computing structural similarities between protein structures using various tools.
    """
    def __init__(self, foldseek_score="evalue"):
        """
        Sets up a tool to compute structural similarities between protein structures.
        """
        self.foldseek_score = foldseek_score
        self.tools = self._set_up_tools()
        self.results = {}    


    def run(self, tool_type: ToolType, pdb_dir: str):
        """
        Run the specified tool on the PDB directory.
        """
        for tool in self.tools:
            if tool.name == tool_type.value:
                print(f"Computing structural similarity with {tool.name}...")
                num_structures = len(os.listdir(pdb_dir))
                # small Gaus without the biggest factor of the sum 
                expected_number_of_scores = (num_structures * (num_structures + 1) // 2) - num_structures
                scores = tool.start_run(pdb_dir)
                if scores is None or len(scores) != expected_number_of_scores:
                    raise ValueError(f"Error: {tool.name} did not return the expected number of scores. Expected {expected_number_of_scores}, got {len(scores) if scores else 'None'}.")
                print(f"Structural similarity computation with {tool.name} completed.")
                return scores
            else:
                continue
        # if no tool matches the tool_type, raise an error
        raise ValueError(f"Tool {tool_type.value} is not available.")


    def _set_up_tools(self):
        """
        Sets up the StructSimTools.
        """
        tools = [foldseek := Foldseek(self.foldseek_score),
                 usalign := USalign(),
                 tmalign := TMalign()
                 ]
        return tools
