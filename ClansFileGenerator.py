import io
from pandas import DataFrame
import pandas as pd
from ClansFile import ClansFile
import os
from utils_for_structures_and_fasta import extract_uids_from_fasta, reset_dir_content
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord



class ClansFileGenerator:
    """
    This class is responsible for generating a CLANS input file.
    It is also able to parse CLANS files and extract the relevant information.
    """
    
    def parse_clans_file(self, clans_file_path: str) -> ClansFile:
        """
        Parses a CLANS file and extracts the relevant information into a ClansFile object.
        Args:
            clans_file_path: A path to the input CLANS file.
        Returns:
            ClansFile: A ClansFile object containing the parsed information.
        """
        print(f"Parsing CLANS file {clans_file_path}...")
        with open(clans_file_path, 'r') as file:
            content = file.read()
        lines = [line.strip() for line in content.strip().splitlines() if line.strip()]
        number_of_sequences = self._parse_number_of_sequences(lines)
        params = self._parse_param_block(lines)
        fasta_records = self._parse_fasta_block(lines)
        coordinates = self._parse_pos_block(lines)
        scores_df = self._parse_scores_block(lines)
        return ClansFile(number_of_sequences, coordinates, scores_df, path_to_fasta=None, fasta_records=fasta_records, parameters=params)


    def _extract_block(self, tag: str, lines: list[str]) -> list[str]|None:
        """
        Extracts a block of text from the lines based on the given tag.
        Args:
            tag (str): The tag indicating the start and end of the block.
            lines (list): The list of lines to search within.

        Returns:
            list: The lines within the extracted block.
        """
        try:
            start = lines.index(f"<{tag}>") + 1
            end = lines.index(f"</{tag}>")
            return lines[start:end]
        except ValueError:
            return None
    
    
    def _parse_number_of_sequences(self, lines: list[str]) -> int:
        """
        Extracts the number of sequences from the lines.
        Args:
            lines (list): The list of lines to search within.
        Returns:
            int: The number of sequences.
        Raises:
            ValueError: If the 'sequences=' line is missing.
        """
        seq_count_line = next((l for l in lines if l.startswith("sequences=")), None)
        if seq_count_line is None:
            raise ValueError("Missing 'sequences=' line.")
        number_of_sequences = int(seq_count_line.split("=")[1])
        return number_of_sequences
    
    
    def _parse_param_block(self, lines: list[str]) -> dict|None:
        """
        Extracts and parses the param block from the lines.
        Args:
            lines (list): The list of lines to search within.
        Returns:
            dict: A dictionary containing the parsed parameters.
        """
        param_block = self._extract_block("param", lines)
        if param_block is None:
            return None
        params = {}
        for line in param_block:
            if not line:
                continue
            key, value = line.split("=")
            params[key.strip()] = value.strip()
        return params
    
    
    def _parse_fasta_block(self, lines: list[str]) -> list[SeqRecord]:
        """
        Extracts and parses the FASTA block from the lines.
        Args:
            lines (list): The list of lines to search within.
        Returns:
            list: A list of SeqRecord objects parsed from the FASTA block.
        """
        seq_block = self._extract_block("seq", lines)
        if seq_block is None:
            raise ValueError("Missing 'seq' block in CLANS file.")
        fasta_text = "\n".join(seq_block)
        fasta_io = io.StringIO(fasta_text)
        fasta_records = list(SeqIO.parse(fasta_io, "fasta"))
        return fasta_records


    def _parse_pos_block(self,lines: list[str]) -> list[tuple[int, float, float, float]]:
        """
        Extracts and parses the position block from the lines.
        Args:
            lines (list): The list of lines to search within.
        Returns:
            list: A list of tuples containing the position information.
        """
        pos_block = self._extract_block("pos", lines)
        if pos_block is None:
            raise ValueError("Missing 'pos' block in CLANS file.")
        positions = []
        for line in pos_block:
            parts = line.split()
            if len(parts) != 4:
                raise ValueError(f"Invalid position line: {line}")
            idx, x, y, z = parts
            positions.append((int(idx), float(x), float(y), float(z)))
        return positions
    
    
    def _parse_scores_block(self, lines: list[str]) -> DataFrame:
        """
        Extracts and parses the scores block from the lines.
        Args:
            lines (list): The list of lines to search within.
        Returns:
            pd.DataFrame: A DataFrame containing the parsed scores.
        """
        hsp_block = self._extract_block("hsp", lines)
        if hsp_block is None:
            raise ValueError("Missing 'hsp' block in CLANS file.")
        scores_data = []
        for line in hsp_block:
            if not line:
                continue
            left, right = line.split(":")
            pdb1, pdb2 = left.split()
            score = right
            scores_data.append({"PDBchain1": int(pdb1), "PDBchain2": int(pdb2), "score": float(score)})
        scores_df = pd.DataFrame(scores_data)
        return scores_df
    

    def generate_clans_file(self, scores: DataFrame, path_to_fasta: str, out_path: str) -> str:
        """
        Generates a CLANS file from a fasta file and the pairwise similarity scores of its sequences.
        Args:
            scores: A pandas DataFrame containing the pairwise similarity scores.
            path_to_fasta: A path to the input fasta file containing the sequences.
            out_path: Path to clans file to be generated.
        Returns:
            The path to the generated CLANS file.
        """
        print(f"Generating CLANS file f{out_path}...")
        uids = extract_uids_from_fasta(path_to_fasta)
        scores = self._transform_scores_to_clans_format(scores, uids)
        self.length_of_fasta = len(uids)
        clans_file = ClansFile(
            self.length_of_fasta,
            self._generate_random_coordinates(self.length_of_fasta),
            scores,
            path_to_fasta=path_to_fasta,
            fasta_records=None,
            parameters=None
        )
        content = clans_file.__str__
        with open(out_path, 'w') as file:
            file.write(content())
        print(f"CLANS file generated at {out_path}")
        return out_path


    def _transform_scores_to_clans_format(self, scores: DataFrame, uids: list[str]) -> DataFrame:
        """
        Transform pairwise similarity scores into the format required by CLANS.
        1. Changes the names of the structures to their corresponding indices of the input fasta.
        2. Changes the order of the rows to this format (with f.e. 3 fasta entries):
         PDBchain1  PDBchain2   score
            0     1      0.8
            0     2      0.7
            1     2      0.6
        3. Drops duplicate lines (f.e. 0 1 0.8 and 1 0 0.8).
        (Duplicate lines should not be in the scores because the scores can still be different for PDBchain1 with PDBchain2
        and PDBchain2 with PDBchain1.)
        
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
            
        Args:
            scores: A pandas DataFrame containing the pairwise similarity scores.
            uids: A list of UIDs in the order they appear in the input fasta file.
        Returns:
            A pandas DataFrame containing the transformed pairwise similarity scores.
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


    def _generate_random_coordinates(self, n) -> list[tuple[int, float, float, float]]:
        """
        Generates n tuples (index, x, y, z) with x, y, z being random floats between [0, 1].
        Args:
            n: Number of tuples to generate.
        Returns:
            A list of tuples (index, x, y, z).
        """
        import random
        tuples = []
        for i in range(n):
            x = random.uniform(0, 1)
            x = f"{x:.3f}"
            y = random.uniform(0, 1)
            y = f"{y:.3f}"
            z = random.uniform(0, 1)
            z = f"{z:.3f}"
            tuples.append((i, float(x), float(y), float(z)))
        return tuples
    