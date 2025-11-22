from StructSimTool import StructSimTool
import os
import pandas as pd
from io import StringIO
from utils_for_PDB import reset_dir_content


class USalign(StructSimTool):
    """
    This class extends the StructSimTool class to implement the USalign tool for protein structure comparison.
    """
    def __init__(self):
        description = "A tool for protein structure comparison using USalign."
        if os.path.exists("USalign_working_dir"):
            reset_dir_content("USalign_working_dir")
        else:
            os.makedirs("USalign_working_dir")
        working_dir = "USalign_working_dir/"
        super().__init__("USalign", description, working_dir)
        self.flag_dir = "-dir" # specifies the directory containing PDB files
        self.flag_outfmt = "-outfmt" # specifies the output format
        self.outfmt = "2"
        
        
    def start_run(self, pdb_dir):
        """
        Initializes the self.command list with the necessary parameters to run the tool and then returns _execute_run with the specified pdb_dir.
        """
        # prepare files for USalign
        pdb_files = os.listdir(pdb_dir)
        pdb_names_list = self.working_dir + "pdb_names.txt"
        with open(pdb_names_list, 'w') as f:
            for pdb_name in pdb_files:
                f.write(f"{pdb_name}\n")
        # example command: "USalign -dir PDBs_for_benchmark/ pdb_names.txt -outfmt 2"
        self.command = [self.name, self.flag_dir, pdb_dir + "/", pdb_names_list, self.flag_outfmt, self.outfmt]
        return self._execute_run()

    
    def _parse_output(self):
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
                # clean PDBchain names
                df2['PDBchain1'] = df2['PDBchain1'].str.split(".").str[0]
                df2['PDBchain2'] = df2['PDBchain2'].str.split(".").str[0]
                return df2
            else:
                print("Failed to parse output: DataFrame is empty.")
                return False
        except Exception as e:
            print(f"Error parsing output: {e}")
            return False
        