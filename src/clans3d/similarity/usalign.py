import logging
import os
from io import StringIO

import pandas as pd

from clans3d.similarity.struct_sim_tool import StructSimTool
from clans3d.utils.file_utils import reset_dir_content

logger = logging.getLogger(__name__)


class USalign(StructSimTool):
    """
    This class extends the StructSimTool class to implement the USalign tool for protein structure comparison.
    """
    def __init__(self, working_dir: str = os.path.join("work", "usalign")):
        description = "A tool for protein structure comparison using USalign."
        super().__init__("USalign", description, working_dir)
        self.flag_dir = "-dir" # specifies the directory containing structure files
        self.flag_outfmt = "-outfmt" # specifies the output format
        self.outfmt_value = "2" # output format 2 is tab-separated with columns: PDBchain1, PDBchain2, TM1, TM2, ...
        
        
    def start_run(self, structures_dir: str, expected_number_of_scores: int) -> pd.DataFrame:
        """
        Initializes the self.command list with the necessary parameters to run the tool and then returns _execute_run with the specified structures_dir.
        Returns the computed similarity scores as a DataFrame.
        
        Args:
            structures_dir (str): The directory containing the structure files to be compared.
            expected_number_of_scores (int): The expected number of scores to be returned by the tool, used for logging purposes.
        Returns:
            pd.DataFrame: A DataFrame containing the similarity scores between the structures.
        """
        # clean working directory before running
        reset_dir_content(self.working_dir)        
        # prepare files for USalign
        structure_files = os.listdir(structures_dir)
        structure_names_list = os.path.join(self.working_dir, "structure_names.txt")
        with open(structure_names_list, 'w') as f:
            for structure_name in structure_files:
                f.write(f"{structure_name}\n")
        # example command: "USalign -dir structures_dir/ structure_names.txt -outfmt 2"
        self.command = [self.name, self.flag_dir, structures_dir + "/", structure_names_list, self.flag_outfmt, self.outfmt_value]
        return self._execute_run(expected_number_of_scores)

    
    def _log_progress(self, stdout_reader: dict,
                      stderr_reader: dict) -> None:
        """Log USalign progress by counting completed alignment pairs.

        USalign writes one TSV data line per pair to stdout (plus a
        header and occasional warnings).  This method counts the data
        lines read so far and logs a percentage based on
        ``self.expected_number_of_scores``.

        Args:
            stdout_reader: Reader dict for the stdout pipe.
            stderr_reader: Reader dict for the stderr pipe.
        """
        pairs_done = self._count_data_lines(stdout_reader)
        if self.expected_number_of_scores > 0:
            percent = min(100, (pairs_done / self.expected_number_of_scores) * 100)
            logger.info(f"{self.name}: Aligning pairs – {int(percent)}% ({pairs_done}/{self.expected_number_of_scores})")
        else:
            logger.info(f"{self.name}: Aligned {pairs_done} pairs so far")


    @staticmethod
    def _count_data_lines(stdout_reader: dict) -> int:
        """Count USalign result lines (excluding header and warnings).

        Args:
            stdout_reader: Reader dict for the stdout pipe.

        Returns:
            Number of data lines read so far.
        """
        count = 0
        for line in stdout_reader["chunks"]:
            stripped = line.strip()
            if stripped and not stripped.startswith("#") and not stripped.startswith("Warning"):
                count += 1
        return count


    def _parse_output(self) -> pd.DataFrame:
        """
        Parses the output of the tool to extract the similarity scores.

        Raises:
            RuntimeError: If the output cannot be parsed or is empty.
        """
        try:
            df = pd.read_csv(StringIO(self.output), sep='\t')
            if df.empty:
                raise RuntimeError(f"{self.name} produced empty output.")
            df.columns = [col.lstrip('#') for col in df.columns]
            df1 = df[['PDBchain1', 'PDBchain2', 'TM1', 'TM2']].copy()
            # take the maximum of TM1 and TM2
            df1['TM'] = df1[['TM1', 'TM2']].max(axis=1)
            # transform TM-score to distance metric
            df1["score"] = 1 - df1["TM"]
            df2 = df1.drop(columns=['TM1', 'TM2'])
            # clean structure names (remove file extensions)
            df2['PDBchain1'] = df2['PDBchain1'].str.split(".").str[0]
            df2['PDBchain2'] = df2['PDBchain2'].str.split(".").str[0]
            # rename columns to generic naming
            df2 = df2.rename(columns={'PDBchain1': 'Sequence_ID_1', 'PDBchain2': 'Sequence_ID_2'})
            return df2
        except RuntimeError:
            raise
        except Exception as e:
            raise RuntimeError(f"{self.name} failed to parse output: {e}") from e
