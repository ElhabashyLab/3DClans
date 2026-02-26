import logging
import os
import subprocess

import pandas as pd

from clans3d.similarity.struct_sim_tool import StructSimTool
from clans3d.utils.file_utils import reset_dir_content

logger = logging.getLogger(__name__)


class Foldseek(StructSimTool):
    """
    This class extends the StructSimTool class to implement the Foldseek tool for protein structure comparison.
    """
    def __init__(self, score):
        description = "A tool for protein structure comparison using Foldseek."
        self.score = score
        working_dir = os.path.join("work", "foldseek")
        super().__init__("foldseek", description, working_dir)
        self.createdb = "createdb"  # command to create a database
        self.search = "search"  # command to search in the database
        self.convertalis = "convertalis"  # command to convert alignment results
        self.flag_format_output = "--format-output" # flag to specify output format
        self.output_columns = self._get_output_columns_from_self_score()  # specifies the columns foldseek returns based on the needed metric (important metrics: qtmscore, ttmscore, evalue)
        self.flag_a = "-a" # enables foldseek to store alignment backtrace information (needed for TM scores)
        self.exhaustive_search = "--exhaustive-search" # enables exhaustive search (slower but more accurate/skips prefiltering)
        self.e = "-e" # e-value threshold for foldseek search
        

    def _get_output_columns_from_self_score(self):
        """
        Generates the ouptut columns of foldseek based on self.score.
        """
        if self.score == "TM":
            output_columns = "query,target,qtmscore,ttmscore"
        elif self.score == "evalue":
            output_columns = "query,target,evalue"
        else:
            raise ValueError(f"Invalid score type: {self.score}. Supported types are 'TM' and 'evalue'.")
        return output_columns


    def start_run(self, structures_dir: str) -> pd.DataFrame:
        """
        Initializes the self.command list with the necessary parameters to run the tool and then returns _execute_run.
        This method should be overridden by subclasses to implement specific tool logic.
        The command structure has 3 steps:
        
        (creating a database)
        foldseek createdb <i:directory|.tsv>|<i:PDB|mmCIF[.gz]|tar[.gz]|DB> ... <i:PDB|mmCIF[.gz]|tar|DB> <o:sequenceDB> [options][/+][-]
        
        (running alignment)
        foldseek search <i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir> [options]
        
        (converting into readable output)
        foldseek convertalis <i:queryDb> <i:targetDb> <i:alignmentDB> <o:alignmentFile> [options][/+][-]
        
        Args:
            structures_dir (str): The directory containing the structure files to be compared.
        Returns:
            pd.DataFrame: A DataFrame containing the similarity scores between the structures.
        """
        # clean working directory before running
        reset_dir_content(self.working_dir)
        # creating target_db = query_db
        self.db = self._create_database(structures_dir, "foldseekDB")
        if not self.db:
            raise Exception("Failed to create Foldseek database.")
        # running alignment
        self.alignmentDb = os.path.join(self.working_dir, "alignmentDb")
        self.tmpDir = os.path.join(self.working_dir, "tmpDir")
        self.command = [self.name, self.search, self.db, self.db, self.alignmentDb, self.tmpDir, self.flag_a] 
        return self._execute_run()


    def _create_database(self, structures_dir: str, db_name: str):
        """
        Creates a foldseek database named db_name from the given structures directory.
        """
        create_db_command = [self.name, self.createdb, structures_dir, os.path.join(self.working_dir, db_name)]
        try:
            result = subprocess.run(create_db_command,
                                    capture_output=True,
                                    text=True,
                                    check=True)
            return os.path.join(self.working_dir, db_name)
        except subprocess.CalledProcessError as e:
            logger.error("Error running %s with %s: %s", self.name, create_db_command, e)
            logger.debug("stdout: %s", e.stdout)
            logger.debug("stderr: %s", e.stderr)
            return False


    def _parse_output(self) -> pd.DataFrame:
        """
        Parses the output of the tool to extract self.score.
        """
        self.alignmentFile = os.path.join(self.working_dir, "alignmentFile")
        convertalis_command = [self.name, self.convertalis, self.db, self.db, self.alignmentDb, self.alignmentFile, self.flag_format_output, self.output_columns]
        try:
            result = subprocess.run(convertalis_command,
                                    capture_output=True,
                                    text=True,
                                    check=True)
            df = pd.read_csv(self.alignmentFile, sep="\t")
            # name the columns of the df
            df.columns = ["Sequence_ID_1", "Sequence_ID_2", "TM1", "TM2"] if self.score == "TM" else ["Sequence_ID_1", "Sequence_ID_2", "evalue"]
            df1 = self._clean_scores(df)
            return df1
        except subprocess.CalledProcessError as e:
            logger.error("Error running %s with %s: %s", self.name, convertalis_command, e)
            logger.debug("stdout: %s", e.stdout)
            logger.debug("stderr: %s", e.stderr)
            return pd.DataFrame()  # Return an empty DataFrame on error


    def _clean_scores(self, df: pd.DataFrame) -> pd.DataFrame:
        """
        Receives a DataFrame containing the columns 'Sequence_ID_1, Sequence_ID_2, score' and returns a cleaned DataFrame with no duplicate rows.
        The score is transformed based on self.score (if self.score is "TM", the maximum of TM1 and TM2 is taken and transformed to a distance metric by 1 - maxTM, if self.score is "evalue", the evalue is taken as score).
        
        Args:
            df (pd.DataFrame): A DataFrame containing the columns 'Sequence_ID_1, Sequence_ID_2, score'
        Returns:
            pd.DataFrame: A cleaned DataFrame with no duplicate rows.
        """
        # remove duplicates like A:B and B:A
        df["pairs"] = df.apply(lambda row: tuple(sorted([row['Sequence_ID_1'], row['Sequence_ID_2']])), axis=1)
        df1 = df.drop_duplicates(subset="pairs")
        df2 = df1.drop(columns=["pairs"])
        # remove rows where A is the same as B
        df3 = df2[df2['Sequence_ID_1'] != df2['Sequence_ID_2']].copy()
        # clean scores based on self.score
        if self.score == "TM":
            df3['maxTM'] = df3[['TM1', 'TM2']].max(axis=1)
            df3["score"] = 1 - df3["maxTM"]
            df4 = df3.drop(columns=['TM1', 'TM2', 'maxTM'])
        else:
            df3['score'] = df3['evalue']
            df4 = df3.drop(columns=['evalue'])
        # reset index
        df4 = df4.reset_index(drop=True)
        return df4
