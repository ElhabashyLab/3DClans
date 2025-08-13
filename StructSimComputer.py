from fasta2PDB import extract_uids_from_fasta, fetch_pdbs_from_uids
from USalign import USalign
from TMalign import TMalign
from Foldseek import Foldseek
from ToolType import ToolType


class StructSimComputer:
    """
    This class is responsible for computing structural similarities between protein structures using various tools.
    """
    def __init__(self):
        """
        Sets up a tool to compute structural similarities between protein structures.
        """
        self.tools = self._set_up_tools()
        self.results = {}    
    
    
    def run(self, tool_type: ToolType, pdb_dir: str):
        """
        Run the specified tool on the PDB directory.
        """
        for tool in self.tools:
            if tool.name == tool_type.value:
                print(f"Computing structural similarity with {tool.name}...")
                scores = tool.start_run(pdb_dir)
                return scores
            else:
                continue
        # if no tool matches the tool_type, raise an error
        raise ValueError(f"Tool {tool_type.value} is not available.")


    def _set_up_tools(self):
        """
        Sets up the StructSimTools.
        """
        tools = [foldseek := Foldseek(),
                 usalign := USalign(),
                 tmalign := TMalign()
                 ]
        return tools
