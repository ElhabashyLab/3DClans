import logging
import os
from pandas import DataFrame
from clans3d.similarity.usalign import USalign
from clans3d.similarity.foldseek import Foldseek
from clans3d.similarity.tool_type import ToolType
from clans3d.similarity.struct_sim_tool import StructSimTool

logger = logging.getLogger(__name__)


class StructSimComputer:
    """
    This class is responsible for computing structural similarities between protein structures using various tools.
    """
    def __init__(self, foldseek_score: str = "evalue", working_dir: str = "work"):
        """
        Sets up a tool to compute structural similarities between protein structures.
        
        Args:
            foldseek_score: Score type for Foldseek ("evalue" or "TM").
            working_dir: Base directory for tool working files. Each tool creates its own subdirectory.
        """
        self.foldseek_score = foldseek_score
        self.working_dir = working_dir
        self.results = {}    


    def run(self, tool_type: ToolType, structures_dir: str) -> DataFrame:
        """
        Run the specified tool on the structures directory.
        """
        tool = self._create_tool(tool_type)
        num_structures = len(os.listdir(structures_dir))
        expected_number_of_scores = (num_structures * (num_structures + 1) // 2) - num_structures
        scores = tool.start_run(structures_dir, expected_number_of_scores)
        number_scores = len(scores)
        if number_scores != expected_number_of_scores:
            logger.warning("%s did not return the expected number of scores. Expected %d, got %d.", tool.name, expected_number_of_scores, number_scores)
        logger.info("Structural similarity computation with %s completed.", tool.name)
        return scores
        

    def _create_tool(self, tool_type: ToolType) -> StructSimTool:
        """Factory method to create the appropriate tool instance."""
        if tool_type == ToolType.FOLDSEEK:
            return Foldseek(self.foldseek_score, os.path.join(self.working_dir, "foldseek"))
        elif tool_type == ToolType.USALIGN:
            return USalign(os.path.join(self.working_dir, "usalign"))
        else:
            raise ValueError(f"Tool {tool_type.value} is not available.")
