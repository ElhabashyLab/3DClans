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
    def __init__(self):
        description = "A tool for protein structure comparison using USalign."
        working_dir = os.path.join("work", "usalign")
        super().__init__("USalign", description, working_dir)
        self.flag_dir = "-dir" # specifies the directory containing structure files
        self.flag_outfmt = "-outfmt" # specifies the output format
        self.outfmt = "2"
        
        
    def start_run(self, structures_dir: str) -> pd.DataFrame:
        """
        Initializes the self.command list with the necessary parameters to run the tool and then returns _execute_run with the specified structures_dir.
        
        Args:
            structures_dir (str): The directory containing the structure files to be compared.
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
        self.command = [self.name, self.flag_dir, structures_dir + "/", structure_names_list, self.flag_outfmt, self.outfmt]
        return self._execute_run()

    
    def _parse_output(self) -> pd.DataFrame:
        """
        Parses the output of the tool to extract the similarity scores.
        """
        try:
            df = pd.read_csv(StringIO(self.output), sep='\t')
            if not df.empty:
                df.columns = [col.lstrip('#') for col in df.columns]
                df1 = df[['PDBchain1', 'PDBchain2', 'TM1', 'TM2']].copy()
                # take the maximum of TM1 and TM2
                df1['TM'] = df1[['TM1', 'TM2']].max(axis=1)
                df1["score"] = 1 - df1["TM"]  # transform TM-score to distance metric
                df2 = df1.drop(columns=['TM1', 'TM2'])
                # clean structure names (remove file extensions)
                df2['PDBchain1'] = df2['PDBchain1'].str.split(".").str[0]
                df2['PDBchain2'] = df2['PDBchain2'].str.split(".").str[0]
                # rename columns to standard naming
                df2 = df2.rename(columns={'PDBchain1': 'Sequence_ID_1', 'PDBchain2': 'Sequence_ID_2'})
                return df2
            else:
                logger.error("Failed to parse output: DataFrame is empty.")
                return pd.DataFrame()  # Return an empty DataFrame if the output is empty
        except Exception as e:
            logger.error("Error parsing output: %s", e)
            return pd.DataFrame()  # Return an empty DataFrame on parsing error
