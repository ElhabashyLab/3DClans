import logging
import os

from clans3d.similarity.usalign import USalign
from clans3d.similarity.tmalign import TMalign
from clans3d.similarity.foldseek import Foldseek
from clans3d.similarity.tool_type import ToolType
from clans3d.similarity.struct_sim_tool import StructSimTool

logger = logging.getLogger(__name__)


class StructSimComputer:
    """
    This class is responsible for computing structural similarities between protein structures using various tools.
    """
    def __init__(self, foldseek_score="evalue"):
        """
        Sets up a tool to compute structural similarities between protein structures.
        """
        self.foldseek_score = foldseek_score
        self.results = {}    


    def run(self, tool_type: ToolType, pdb_dir: str):
        """
        Run the specified tool on the PDB directory.
        """
        tool = self._create_tool(tool_type)
        logger.info("Computing structural similarity with %s...", tool.name)
        num_structures = len(os.listdir(pdb_dir))
        expected_number_of_scores = (num_structures * (num_structures + 1) // 2) - num_structures
        scores = tool.start_run(pdb_dir)
        if scores is None or len(scores) != expected_number_of_scores:
            if scores is None:
                len_scores = 0
            else:
                len_scores = len(scores)
            logger.warning("%s did not return the expected number of scores. Expected %d, got %d.", tool.name, expected_number_of_scores, len_scores)
        logger.info("Structural similarity computation with %s completed.", tool.name)
        return scores
        

    def _create_tool(self, tool_type: ToolType) -> StructSimTool:
        """Factory method to create the appropriate tool instance."""
        if tool_type == ToolType.FOLDSEEK:
            return Foldseek(self.foldseek_score)
        elif tool_type == ToolType.USALIGN:
            return USalign()
        elif tool_type == ToolType.TMALIGN:
            return TMalign()
        else:
            raise ValueError(f"Tool {tool_type.value} is not available.")
