import io
from pandas import DataFrame
import pandas as pd
from clans3d.core.clans_file import ClansFile
from clans3d.utils.fasta_utils import extract_uids_from_fasta, extract_uid_from_recordID
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
        The indices in the CLANS file are mapped back to the actual sequence UIDs based on
        the order in the <seq> block.
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
        coordinates_raw = self._parse_pos_block(lines)
        scores_df_raw = self._parse_scores_block(lines)
        
        # Build index-to-UID mapping based on fasta_records order
        # Use short UID (via extract_uid_from_recordID) to match format used throughout the system
        idx_to_uid = {
            i: extract_uid_from_recordID(record.id) if record.id else f"seq_{i}" 
            for i, record in enumerate(fasta_records)
        }
        
        # Map indices back to UIDs in coordinates
        coordinates: list[tuple[str, float, float, float]] = [
            (idx_to_uid[idx], x, y, z)
            for idx, x, y, z in coordinates_raw
        ]
        
        # Map indices back to UIDs in scores DataFrame
        scores_df = scores_df_raw.copy()
        scores_df["PDBchain1"] = scores_df["PDBchain1"].map(idx_to_uid)
        scores_df["PDBchain2"] = scores_df["PDBchain2"].map(idx_to_uid)
        
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
        Scores should contain UIDs (not indices) - the conversion to indices happens in ClansFile.__str__.
        Args:
            scores: A pandas DataFrame containing the pairwise similarity scores with UIDs.
            path_to_fasta: A path to the input fasta file containing the sequences.
            out_path: Path to clans file to be generated.
        Returns:
            The path to the generated CLANS file.
        """
        print(f"Generating CLANS file {out_path}...")
        uids = extract_uids_from_fasta(path_to_fasta)
        scores = self._normalize_scores_format(scores, uids)
        self.length_of_fasta = len(uids)
        clans_file = ClansFile(
            self.length_of_fasta,
            self._generate_random_coordinates(uids),
            scores,
            path_to_fasta=path_to_fasta,
            fasta_records=None,
            parameters=None
        )
        content = str(clans_file)
        with open(out_path, 'w') as file:
            file.write(content)
        print(f"CLANS file generated at {out_path}")
        return out_path


    def _normalize_scores_format(self, scores: DataFrame, uids: list[str]) -> DataFrame:
        """
        Normalize pairwise similarity scores format:
        1. Ensures consistent ordering of pairs (uid1 < uid2 alphabetically for deduplication)
        2. Drops duplicate pairs
        3. Keeps UIDs as identifiers (conversion to indices happens in ClansFile.__str__)
        
        Args:
            scores: A pandas DataFrame containing the pairwise similarity scores with columns
                    [PDBchain1, PDBchain2, score] where PDBchain1/2 are UIDs.
            uids: A list of UIDs in the order they appear in the input fasta file.
        Returns:
            A pandas DataFrame containing the normalized pairwise similarity scores with UIDs.
        """
        scores1 = scores.copy()
        # Ensure consistent ordering for deduplication (smaller UID first alphabetically)
        mask = scores1["PDBchain1"] > scores1["PDBchain2"]
        scores1.loc[mask, ["PDBchain1", "PDBchain2"]] = scores1.loc[mask, ["PDBchain2", "PDBchain1"]].values
        # Drop duplicates
        scores1 = scores1.drop_duplicates(subset=["PDBchain1", "PDBchain2"]).reset_index(drop=True)
        return scores1


    def _generate_random_coordinates(self, uids: list[str]) -> list[tuple[str, float, float, float]]:
        """
        Generates tuples (uid, x, y, z) with x, y, z being random floats between [0, 1].
        Args:
            uids: List of sequence UIDs.
        Returns:
            A list of tuples (uid, x, y, z).
        """
        import random
        tuples = []
        for uid in uids:
            x = float(f"{random.uniform(0, 1):.3f}")
            y = float(f"{random.uniform(0, 1):.3f}")
            z = float(f"{random.uniform(0, 1):.3f}")
            tuples.append((uid, x, y, z))
        return tuples
