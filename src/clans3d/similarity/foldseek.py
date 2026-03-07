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
    
    How to run Foldseek:
    The command structure has 3 steps ->
        
        1. (create a database)
        foldseek createdb <i:directory|.tsv>|<i:PDB|mmCIF[.gz]|tar[.gz]|DB> ... <i:PDB|mmCIF[.gz]|tar|DB> <o:sequenceDB> [options][/+][-]
        
        2. (run alignment)
        foldseek search <i:queryDB> <i:targetDB> <o:alignmentDB> <tmpDir> [options]
        
        3. (convert into readable output)
        foldseek convertalis <i:queryDb> <i:targetDb> <i:alignmentDB> <o:alignmentFile> [options][/+][-]
    """
    def __init__(self, score):
        description = "A tool for protein structure comparison using Foldseek."
        self.score = score
        working_dir = os.path.join("work", "foldseek")
        super().__init__("foldseek", description, working_dir)
        self.command_createdb = "createdb"  # command to create a database
        self.command_search = "search"  # command to search in the database
        self.command_convertalis = "convertalis"  # command to convert alignment results
        self.flag_format_output = "--format-output" # flag to specify output format
        self.arg_output_columns = self._get_output_columns_from_self_score()  # specifies the columns foldseek returns based on the needed metric (important metrics: qtmscore, ttmscore, evalue)
        self.flag_a = "-a" # enables foldseek to store alignment backtrace information (needed for TM scores)
        self.flag_exhaustive_search = "--exhaustive-search" # enables exhaustive search (slower but more accurate/skips prefiltering)
        self.flag_e_value_threshold = "-e" # e-value threshold for foldseek search
        self._last_logged_phase: str = "" # tracks the last logged phase for progress logging


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


    def start_run(self, structures_dir: str, expected_number_of_scores: int) -> pd.DataFrame:
        """
        Initializes a foldseek database with the structures_dir
        and sets self.command list with the necessary parameters to run the tool and then returns _execute_run.
        
        Args:
            structures_dir (str): The directory containing the structure files to be compared.
            expected_number_of_scores (int): The expected number of scores to be returned by the tool, used for logging purposes.
        Returns:
            pd.DataFrame: A DataFrame containing the similarity scores between the structures.
        """
        # clean working directory before running
        reset_dir_content(self.working_dir)
        # creating target_db = query_db
        self.db = self._create_database(structures_dir, "foldseekDB")
        # running alignment
        self.alignmentDb = os.path.join(self.working_dir, "alignmentDb")
        self.tmpDir = os.path.join(self.working_dir, "tmpDir")
        self.command = [self.name, self.command_search, self.db, self.db, self.alignmentDb, self.tmpDir, self.flag_a] 
        return self._execute_run(expected_number_of_scores)


    # ---- Foldseek phase keywords (in execution order) ----
    _PHASE_KEYWORDS: list[tuple[str, str]] = [
        ("counting k-mers",     "Indexing (counting k-mers) ..."),
        ("index table: fill",   "Indexing (fill) ..."),
        ("prefiltering",        "Prefiltering ..."),
        ("structurealign",      "Structure alignment ..."),
        ("convertalis",         "Converting results ..."),
    ]

    def _log_progress(self, stdout_reader: dict,
                      stderr_reader: dict) -> None:
        """Log Foldseek progress by detecting phase transitions.

        When piped, Foldseek writes all output (including phase labels)
        to stdout and nothing to stderr.  Progress bars only appear as
        a final "done" line per phase, so instead of percentages this
        method detects phase transitions and logs them once.

        Args:
            stdout_reader: Reader dict for the stdout pipe.
            stderr_reader: Reader dict for the stderr pipe.
        """
        phase = self._detect_phase(stdout_reader)

        if phase and phase != self._last_logged_phase:
            self._last_logged_phase = phase
            logger.info(f"{self.name}: {phase}")


    def _detect_phase(self, reader: dict) -> str:
        """Determine the current Foldseek processing phase.

        Scans the captured stdout lines (from newest to oldest) for
        known phase keywords and returns a human-readable label.

        Args:
            reader: Reader dict for the stdout pipe.

        Returns:
            A descriptive phase label, e.g. ``"Prefiltering"``,
            or an empty string if no phase has been detected yet.
        """
        for raw_line in reversed(reader["chunks"]):
            line_lower = raw_line.lower()
            for keyword, label in self._PHASE_KEYWORDS:
                if keyword in line_lower:
                    return label
        return ""


    def _create_database(self, structures_dir: str, db_name: str) -> str:
        """
        Creates a foldseek database named db_name from the given structures directory.
        Args:
            structures_dir (str): The directory containing the structure files to be compared.
            db_name (str): The name of the database to be created.
        Returns:
            The path to the created database, or False if the database creation failed.
        """
        create_db_command = [self.name, self.command_createdb, structures_dir, os.path.join(self.working_dir, db_name)]
        try:
            subprocess.run(create_db_command,
                           capture_output=True,
                           text=True,
                           check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"foldseek createdb failed (exit {e.returncode}): {e.stderr}"
            ) from e
        return os.path.join(self.working_dir, db_name)


    def _parse_output(self) -> pd.DataFrame:
        """
        Parses the output of the tool to extract self.score.
        """
        self.alignmentFile = os.path.join(self.working_dir, "alignmentFile")
        convertalis_command = [self.name, self.command_convertalis, self.db, self.db, self.alignmentDb, self.alignmentFile, self.flag_format_output, self.arg_output_columns]
        try:
            subprocess.run(convertalis_command,
                           capture_output=True,
                           text=True,
                           check=True)
        except subprocess.CalledProcessError as e:
            raise RuntimeError(
                f"foldseek convertalis failed (exit {e.returncode}): {e.stderr}"
            ) from e
        df = pd.read_csv(self.alignmentFile, sep="\t")
        if df.empty:
            raise RuntimeError(f"{self.name} convertalis produced an empty output file.")
        # rename the columns of the df
        df.columns = ["Sequence_ID_1", "Sequence_ID_2", "TM1", "TM2"] if self.score == "TM" else ["Sequence_ID_1", "Sequence_ID_2", "evalue"]
        return self._clean_scores(df)


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
