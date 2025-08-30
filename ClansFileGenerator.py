from ClansFile import ClansFile
import os
from fasta2PDB import extract_uids_from_fasta, delete_dir_content


class ClansFileGenerator:
    """
    This class is responsible for generating a CLANS input file from a fasta file.
    """
    def __init__(self):
        self.output_dir = "clans_files"
        if os.path.exists(self.output_dir):
            delete_dir_content(self.output_dir)
        else:
            os.makedirs(self.output_dir)


    def generate_clans_file(self, scores, fasta):
        """
        Generates a CLANS input file from a fasta file and the pairwise similarity scores of the sequences.
        Returns the path to the generated CLANS file.
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
        clans_file_path = os.path.join(self.output_dir, "file.clans")
        with open(clans_file_path, 'w') as file:
            file.write(content())
        print(f"CLANS file generated at {clans_file_path}")
        return clans_file_path


    def _transform_scores_to_clans_format(self, scores, uids):
        """
        Transform the pairwise similarity scores into the format required by CLANS.
        1. Changes the names of the PDBs to their corresponding indices of the input fasta.
        2. Transforms the TM-score to a distance metric (1 - TM-score)
        3. Changes the order of the rows to this format (with f.e. 3 fasta entries):
         PDBchain1  PDBchain2   TM
            0     1      0.8
            0     2      0.7
            1     2      0.6
        
        example:
        
        order of sequences in fasta file:
            0: A0A2M8UWB6
            1: A0A836ZK00
            2: A0A2W5V1G2
        scores:
            PDBchain1  PDBchain2   TM
            A0A836ZK00 A0A2W5V1G2 0.4469
            A0A836ZK00 A0A2M8UWB6 0.4481
            A0A2M8UWB6 A0A2W5V1G2 0.8499
        transformed scores:
            PDBchain1  PDBchain2   TM
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
        # transform TM-score to distance metric
        scores1["TM"] = 1 - scores1["TM"]
        # order the pairs
        scores1 = scores1.sort_values(by=["PDBchain1", "PDBchain2"]).reset_index(drop=True)
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
