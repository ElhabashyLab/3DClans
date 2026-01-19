import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class ClansFile:
    """
    This Class represents a CLANS file
    A Clans file consists of:
    - number of sequences
    - metadata parameters
    - the fasta of the sequences
    - coordinates of the sequences
    - pairwise similarity scores of the sequences
    
    (example):
    sequences=3
    <param>
    pval=1e-5
    </param>
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
    def __init__(self,
                 number_of_sequences: int,
                 coordinates: list[tuple[int, float, float, float]],
                 scores: pd.DataFrame,
                 path_to_fasta: str|None = None,
                 fasta_records: None|list[SeqRecord] = None,
                 parameters: dict|None = None):
        self.number_of_sequences = number_of_sequences
        self.parameters = parameters
        self.path_to_fasta = path_to_fasta
        if fasta_records is not None:
            self.fasta_records = fasta_records
        elif path_to_fasta is not None:
            self.fasta_records = SeqIO.parse(path_to_fasta, "fasta")
        else:
            raise ValueError("Either fasta_records or path_to_fasta must be provided to generate a ClansFile object.")
        self.coordinates = coordinates
        self.scores = scores
        
        
    def __str__(self) -> str:
        """
        Returns a string representation of the CLANS file.
        """
        content = []
        content.append(f"sequences={self.number_of_sequences}")
        content = self._add_parameters_to_content(content)
        content = self._add_fasta_to_content(content)
        content = self._add_coordinates_to_content(content)
        content = self._add_scores_to_content(content)
        return "\n".join(content)


    def _add_parameters_to_content(self, content):
        content.append("<param>")
        if self.parameters is not None:
            for key, value in self.parameters.items():
                content.append(f"-{key} {value}")
        content.append("</param>")
        return content


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
    