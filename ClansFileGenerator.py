from ClansFile import ClansFile
import os
from fasta2PDB import extract_uids_from_fasta, delete_dir_content
import numpy as np


class ClansFileGenerator:
    """
    This class is responsible for generating a CLANS input file from a fasta file.
    It is also able to parse CLANS files and extract the relevant information.
    """
    def __init__(self, output_dir="clans_files"):
        self.output_dir = output_dir
        if os.path.exists(self.output_dir):
            delete_dir_content(self.output_dir)
        else:
            os.makedirs(self.output_dir)

    
    @staticmethod
    def parse_clans_file(clans_file_path):
        """
        Parses a CLANS input file and extracts the relevant information.
        Returns a ClansFile object.
        Args:
            clans_file_path: A path to the input CLANS file.
        """
        print(f"Parsing CLANS file {clans_file_path}...")
        print("debug 3")
        with open(clans_file_path, 'r') as file:
            content = file.read()
        clans_file = ClansFile.from_string(content)
        return clans_file
    

    def generate_clans_file(self, scores, fasta):
        """
        Generates a CLANS input file from a fasta file and the pairwise similarity scores of the sequences.
        Returns the path to the generated CLANS file.
        Args:
            scores: A pandas DataFrame containing the pairwise similarity scores.
            fasta: A path to the input fasta file.
        """
        print(f"Generating CLANS file in {self.output_dir}...")
        uids = extract_uids_from_fasta(fasta)
        scores = self._transform_scores_to_clans_format(scores, uids)
        self.length_of_fasta = len(uids)
        clans_file = ClansFile(
            self.length_of_fasta,
            fasta,
            self._generate_random_coordinates(self.length_of_fasta),
            scores
        )
        content = clans_file.__str__
        fasta_name = os.path.basename(fasta).split(".")[0]
        clans_file_path = os.path.join(self.output_dir, f"{fasta_name}.clans")
        with open(clans_file_path, 'w') as file:
            file.write(content())
        print(f"CLANS file generated at {clans_file_path}")
        return clans_file_path


    def _transform_scores_to_clans_format(self, scores, uids):
        """
        Transform the pairwise similarity scores into the format required by CLANS.
        1. Changes the names of the PDBs to their corresponding indices of the input fasta.
        2. Changes the order of the rows to this format (with f.e. 3 fasta entries):
         PDBchain1  PDBchain2   score
            0     1      0.8
            0     2      0.7
            1     2      0.6
        3. Drops duplicate lines (f.e. 0 1 0.8 and 1 0 0.8).
        
        (Duplicate lines should not be in the scores becauese the scores can still be different for PDBchain1 with PDBchain2
        and PDBchain2 with PDBchain1!))
        
        example:
        
        order of sequences in fasta file:
            0: A0A2M8UWB6
            1: A0A836ZK00
            2: A0A2W5V1G2
        scores:
            PDBchain1  PDBchain2   score
            A0A836ZK00 A0A2W5V1G2 0.4469
            A0A2W5V1G2 A0A836ZK00 0.4469
            A0A836ZK00 A0A2M8UWB6 0.4481
            A0A2M8UWB6 A0A2W5V1G2 0.8499
        transformed scores:
            PDBchain1  PDBchain2   score
            0   1   0.4469
            0   2   0.4481
            1   2   0.8499
        """
        # map uid to index
        uid_to_index = {uid: i for i, uid in enumerate(uids)}
        # replace uids with index in scores
        scores1 = scores.copy()
        scores1["PDBchain1"] = scores1["PDBchain1"].map(lambda x: uid_to_index[x])
        scores1["PDBchain2"] = scores1["PDBchain2"].map(lambda x: uid_to_index[x])
        # make sure PDBchain1 < PDBchain2
        scores1["PDBchain1_new"] = np.minimum(scores1["PDBchain1"], scores1["PDBchain2"])
        scores1["PDBchain2_new"] = np.maximum(scores1["PDBchain1"], scores1["PDBchain2"])
        scores1[["PDBchain1", "PDBchain2"]] = scores1[["PDBchain1_new", "PDBchain2_new"]]
        scores1 = scores1.drop(columns=["PDBchain1_new", "PDBchain2_new"])
        # order the pairs
        scores1 = scores1.sort_values(by=["PDBchain1", "PDBchain2"]).reset_index(drop=True)
        # drop duplicates
        scores1 = scores1.drop_duplicates(subset=["PDBchain1", "PDBchain2"]).reset_index(drop=True)
        return scores1


    def _generate_random_coordinates(self, number_of_sequences):
        """
        Generates random coordinates for the sequences.
        """
        import random
        coordinates = []
        for i in range(number_of_sequences):
            x = random.uniform(0, 1)
            x = f"{x:.3f}"
            y = random.uniform(0, 1)
            y = f"{y:.3f}"
            z = random.uniform(0, 1)
            z = f"{z:.3f}"
            coordinates.append((i, x, y, z))
        return coordinates
