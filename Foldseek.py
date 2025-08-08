from StructSimTool import StructSimTool
import subprocess
import os
import pandas as pd
"""
This class extends the StructSimTool class to implement the Foldseek tool for protein structure comparison.
"""


class Foldseek(StructSimTool):
    def __init__(self):
        description = "A tool for protein structure comparison using Foldseek."
        os.makedirs("Foldseek_working_dir", exist_ok=True)
        working_dir = "Foldseek_working_dir/"
        super().__init__("foldseek", description, working_dir)
        self.createdb = "createdb"  # command to create a database
        self.search = "search"  # command to search in the database
        self.convertalis = "convertalis"  # command to convert alignment results
        self.flag_format_output = "--format-output" # flag to specify output format
        self.output_columns = "query,target,qtmscore,ttmscore"  # columns to output
        self.flag_a = "-a" # enables foldseek to store alignment backtrace information (needed for TM scores)


    def start_run(self, pdb_dir):
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
        """
        # creating target_db = query_db
        self.db = self._create_database(pdb_dir, "foldseekDB")
        if not self.db:
            raise Exception("Failed to create Foldseek database.")
        # running alignment
        self.alignmentDb = self.working_dir + "alignmentDb"
        self.tmpDir = self.working_dir + "tmpDir"
        self.command = [self.name, self.search, self.db, self.db, self.alignmentDb, self.tmpDir, self.flag_a] 
        return self._execute_run()


    def _create_database(self, pdb_dir, db_name):
        """
        Creates a foldseek database named db_name with pdb_dir.
        """
        create_db_command = [self.name, self.createdb, pdb_dir, self.working_dir + db_name]
        try:
            result = subprocess.run(create_db_command,
                                    capture_output=True,
                                    text=True,
                                    check=True)
            return self.working_dir + db_name
        except subprocess.CalledProcessError as e:
            print(f"Error running {self.name} with {create_db_command}: {e}")
            print(f"stdout: {e.stdout}")
            print(f"stderr: {e.stderr}")
            return False


    def _parse_output(self):
        """
        Parses the output of the tool to extract the similarity scores.
        """
        self.alignmentFile = self.working_dir + "alignmentFile"
        convertalis_command = [self.name, self.convertalis, self.db, self.db, self.alignmentDb, self.alignmentFile, self.flag_format_output, self.output_columns]
        try:
            result = subprocess.run(convertalis_command,
                                    capture_output=True,
                                    text=True,
                                    check=True)
            df = pd.read_csv(self.alignmentFile, sep="\t")
            df.columns = ['PDBchain1', 'PDBchain2', 'TM1', 'TM2']
            df = self._clean_scores(df)
            print(df)
            return df
        except subprocess.CalledProcessError as e:
            print(f"Error running {self.name} with {convertalis_command}: {e}")
            print(f"stdout: {e.stdout}")
            print(f"stderr: {e.stderr}")
            return False

    def _clean_scores(self, df):
        """
        recieves a DataFrame containing PDBchain1, PDBchain2, TM1, TM2 columns
        and returns a cleaned DataFrame with no duplicate rows like PDBchain1:PDBchain2 and PDBchain2:PDBchain1
        """
        # remove duplicates like PDBchain1:PDBchain2 and PDBchain2:PDBchain1
        df["pairs"] = df.apply(lambda row: tuple(sorted([row['PDBchain1'], row['PDBchain2']])), axis=1)
        df = df.drop_duplicates(subset="pairs")
        df = df.drop(columns=["pairs"])
        # remove rows where PDBchain1 is the same as PDBchain2
        df = df[df['PDBchain1'] != df['PDBchain2']]
        df['TM'] = df[['TM1', 'TM2']].max(axis=1)
        # reset index
        df = df.reset_index(drop=True)
        return df