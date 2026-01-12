import io
import pandas as pd
from Bio import SeqIO


class ClansFile:
    """
    This Class represents a CLANS file
    A Clans file consists of:
    - number of sequences
    - the fasta of the sequences (can be None if fasta_records are provided)
    - coordinates of the sequences (at first randomly generated)
    - pairwise similarity scores of the sequences
    - the user can alternatively provide the fasta records directly
    
    (example):
    sequences=3
    <seq>
    >tr|A0A2M8UWB6|A0A2M8UWB6_PSESP/524-595
    SPQEFASQLAKPLGAKTAQREVEQALRDLHLPFDERRPYALRRLRDRIEANLSGLMGPSVAQDMVESFLPYK
    >tr|A0A836ZK00|A0A836ZK00_PSESY/73-144
    SPQEFATQLAKPLGAKAAQKEVEQALRDLYLPFDERRPYALRRLRDRIEANLSGLMGPSVAQDMVETFLPYK
    >tr|A0A2W5V1G2|A0A2W5V1G2_PSEOL/284-355
    SPQDFAAQLAKPLGAKTAQREVEQALRDLHLPFDESRPYALRRLRDRIEANLSGLMGPSVAQDIVETFLPYK
    </seq>
    <pos>
    0 0.676 0.371 0.721
    1 0.655 0.385 0.015
    2 0.387 0.890 0.454
    </pos>
    <hsp>
    0 1:4.6e-50
    0 2:4.8e-50
    1 2:2.8e-48
    </hsp>

    """
    def __init__(self, number_of_sequences, path_to_fasta, coordinates, scores, fasta_records=None):
        self.number_of_sequences = number_of_sequences
        self.path_to_fasta = path_to_fasta
        if fasta_records is not None:
            self.fasta_records = fasta_records
        else:
            self.fasta_records = SeqIO.parse(path_to_fasta, "fasta")
        self.coordinates = coordinates
        self.scores = scores
        
        
    def __str__(self):
        """
        Returns a string representation of the CLANS file.
        """
        content = []
        content.append(f"sequences={self.number_of_sequences}")
        content = self._add_fasta_to_content(content)
        content = self._add_coordinates_to_content(content)
        content = self._add_scores_to_content(content)
        return "\n".join(content)
    
    
    def get_coordinates(self):
        """
        Returns a pandas DataFrame with the coordinates of each PDBchain in the Clans file.
        """
        df_coords = pd.DataFrame(self.coordinates, columns=["PDBchain1", "x", "y", "z"])
        df_coords[["x", "y", "z"]] = df_coords[["x", "y", "z"]].apply(pd.to_numeric, errors="coerce")
        df_coords = df_coords.apply(pd.to_numeric, errors='coerce')
        return df_coords


    @classmethod
    def from_string(cls, string: str):
        """
        Parses a CLANS file from a string and returns a ClansFile object.
        Args:
            string: A string representation of a CLANS file.
        Returns:
            ClansFile: A ClansFile object.
        """
        lines = [line.strip() for line in string.strip().splitlines() if line.strip()]
        number_of_sequences = cls._parse_number_of_sequences(lines)
        params = cls._parse_param_block(lines)
        fasta_records = cls._parse_fasta_block(lines)
        coordinates = cls._parse_pos_block(lines)
        scores_df = cls._parse_scores_block(lines)

        return cls(number_of_sequences, None, coordinates, scores_df, fasta_records=fasta_records)
    
    
    @staticmethod
    def _parse_scores_block(lines):
        """Extracts and parses the scores block from the lines.

        Args:
            lines (list): The list of lines to search within.
        Returns:
            pd.DataFrame: A DataFrame containing the parsed scores.
        """
        hsp_block = ClansFile._extract_block("hsp", lines)
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
    
    
    @staticmethod
    def _parse_fasta_block(lines):
        """Extracts and parses the FASTA block from the lines.

        Args:
            lines (list): The list of lines to search within.
        Returns:
            list: A list of SeqRecord objects parsed from the FASTA block.
        """
        seq_block = ClansFile._extract_block("seq", lines)
        fasta_text = "\n".join(seq_block)
        fasta_io = io.StringIO(fasta_text)
        fasta_records = list(SeqIO.parse(fasta_io, "fasta"))
        return fasta_records


    @staticmethod
    def _parse_pos_block(lines):
        """Extracts and parses the position block from the lines.

        Args:
            lines (list): The list of lines to search within.
        Returns:
            list: A list of tuples containing the position information.
        """
        pos_block = ClansFile._extract_block("pos", lines)
        positions = []
        for line in pos_block:
            parts = line.split()
            if len(parts) != 4:
                raise ValueError(f"Invalid position line: {line}")
            idx, x, y, z = parts
            positions.append((int(idx), float(x), float(y), float(z)))
        return positions


    @staticmethod
    def _parse_number_of_sequences(lines):
        """Extracts the number of sequences from the lines.

        Args:
            lines (list): The list of lines to search within.

        Raises:
            ValueError: If the 'sequences=' line is missing.
        """
        seq_count_line = next((l for l in lines if l.startswith("sequences=")), None)
        if seq_count_line is None:
            raise ValueError("Missing 'sequences=' line.")
        number_of_sequences = int(seq_count_line.split("=")[1])
        return number_of_sequences
    
    
    @staticmethod
    def _parse_param_block(lines):
        """Extracts and parses the param block from the lines.

        Args:
            lines (list): The list of lines to search within.
        Returns:
            dict: A dictionary containing the parsed parameters.
        """
        param_block = ClansFile._extract_block("param", lines)
        params = {}
        for line in param_block:
            if not line:
                continue
            key, value = line.split("=")
            params[key.strip()] = value.strip()
        return params


    @staticmethod
    def _extract_block(tag, lines):
        """Extracts a block of text from the lines based on the given tag.

        Args:
            tag (str): The tag indicating the start and end of the block.
            lines (list): The list of lines to search within.

        Returns:
            list: The lines within the extracted block.
        """
        start = lines.index(f"<{tag}>") + 1
        end = lines.index(f"</{tag}>")
        return lines[start:end]
        
    
    def _add_fasta_to_content(self, content):
        content.append("<seq>")
        for record in self.fasta_records:
            header = f">{record.description}"
            sequence = str(record.seq)
            content.append(f"{header}\n{sequence}")
        content.append("</seq>")
        return content
    
    
    def _add_coordinates_to_content(self, content):
        content.append("<pos>")
        for coord in self.coordinates:
            content.append(f"{coord[0]} {coord[1]} {coord[2]} {coord[3]}")
        content.append("</pos>")
        return content
    
    
    def _add_scores_to_content(self, content):
        content.append("<hsp>")
        for index, row in self.scores.iterrows():
            content.append(f"{int(row['PDBchain1'])} {int(row['PDBchain2'])}:{row['score']}")
        content.append("</hsp>")
        return content
    