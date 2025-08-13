from ClansFile import ClansFile
import os
from Bio import SeqIO


class ClansFileGenerator:
    """
    This class is responsible for generating a CLANS input file from a fasta file.
    """
    def __init__(self):
        self.output_dir = "clans_files"
        os.makedirs(self.output_dir, exist_ok=True)


    def generate_clans_file(self, scores, fasta):
        """
        Generates a CLANS input file from a fasta file and the pairwise similarity scores of the sequences.
        Returns the path to the generated CLANS file.
        """
        # need length of fasta file
        self.length_of_fasta = sum(1 for _ in SeqIO.parse(fasta, "fasta"))
        # need content of fasta file
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
        return clans_file_path


    def _generate_random_coordinates(self, number_of_sequences):
        """
        Generates random coordinates for the sequences.
        """
        import random
        coordinates = []
        for i in range(number_of_sequences):
            x = random.uniform(0, 1)
            y = random.uniform(0, 1)
            z = random.uniform(0, 1)
            coordinates.append((i, x, y, z))
        return coordinates
