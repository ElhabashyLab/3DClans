import logging
import subprocess
from clans3d.utils.file_utils import reset_dir_content

logger = logging.getLogger(__name__)


class StructSimTool():
    """
    This module defines a base class for tools that can be used in a benchmarking context.
    """
    def __init__(self, name: str, description: str, working_dir: str):
        self.name = name
        self.description = description
        self.working_dir = working_dir
        self.command = []
        self.output = None
        
        
    def start_run(self, structures_dir):
        """
        Initializes the self.command list with the necessary parameters to run the tool and then returns _execute_run.
        This method should be overridden by subclasses to implement specific tool logic.
        """
        raise NotImplementedError("Subclasses should implement this method to run the tool.")
    
    
    def _execute_run(self):
        """
        Executes self.command in a subprocess and captures the output.
        """
        logger.debug("Running command: %s", " ".join(self.command))
        try:
            result = subprocess.run(self.command,
                                    capture_output=True,
                                    text=True,
                                    check=True)
            self.output = result.stdout
            scores = self._parse_output()
        except subprocess.CalledProcessError as e:
            logger.error("Error running %s with %s: %s", self.name, self.command, e)
            scores = None
        return scores
    
    
    def _parse_output(self):
        """
        Parses the output of the tool to extract the similarity scores.
        This method should be overridden by subclasses to implement specific parsing logic.
        """
        raise NotImplementedError("Subclasses should implement this method to parse output.")
