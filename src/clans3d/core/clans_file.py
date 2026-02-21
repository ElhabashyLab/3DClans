import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from clans3d.utils.fasta_utils import extract_uid_from_recordID


class ClansFile:
    """
    This Class represents a CLANS file
    A Clans file consists of:
    - number of sequences
    - metadata parameters
    - the fasta of the sequences
    - coordinates of the sequences (stored with UIDs, not indices)
    - pairwise similarity scores of the sequences (stored with UIDs, not indices)
    
    (example CLANS file format - note: internally UIDs are used instead of indices):
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
                 coordinates: list[tuple[str, float, float, float]],
                 scores: pd.DataFrame,
                 path_to_fasta: str|None = None,
                 fasta_records: None|list[SeqRecord] = None,
                 parameters: dict|None = None):
        self.number_of_sequences = number_of_sequences
        self.parameters = parameters
        self.path_to_fasta = path_to_fasta
        if fasta_records is not None:
            self.fasta_records = list(fasta_records)
        elif path_to_fasta is not None:
            self.fasta_records = list(SeqIO.parse(path_to_fasta, "fasta"))
        else:
            raise ValueError("Either fasta_records or path_to_fasta must be provided to generate a ClansFile object.")
        self.coordinates = coordinates
        self.scores = scores
        # Build UID to index mapping based on fasta_records order (needed for __str__)
        # Use short UID (extracted via extract_uid_from_recordID) to match the format used in scores/coordinates
        self._uid_to_index = {
            extract_uid_from_recordID(record.id): i 
            for i, record in enumerate(self.fasta_records)
        }
        
        # Validate data consistency
        self._validate_data_consistency()
        
        
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

    
    def _validate_data_consistency(self):
        """
        Validates that the data is consistent:
        - Number of coordinates matches number of sequences
        - All UIDs in coordinates exist in the fasta_records
        - All UIDs in scores exist in the fasta_records
        """
        # Validate coordinates count
        if len(self.coordinates) != self.number_of_sequences:
            raise ValueError(
                f"Number of coordinates ({len(self.coordinates)}) does not match "
                f"number of sequences ({self.number_of_sequences})."
            )
        
        # Validate coordinates UIDs
        coord_uids = {coord[0] for coord in self.coordinates}
        missing_coord_uids = coord_uids - set(self._uid_to_index.keys())
        if missing_coord_uids:
            raise ValueError(
                f"Coordinates contain UIDs not found in fasta_records: {missing_coord_uids}"
            )
        
        # Validate scores UIDs
        if not self.scores.empty:
            score_uids = set(self.scores['PDBchain1'].unique()) | set(self.scores['PDBchain2'].unique())
            missing_score_uids = score_uids - set(self._uid_to_index.keys())
            if missing_score_uids:
                raise ValueError(
                    f"Scores contain UIDs not found in fasta_records: {missing_score_uids}"
                )


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
            # Convert UID back to index for CLANS file format
            uid = coord[0]
            idx = self._uid_to_index.get(uid, uid)  # Fallback to uid if not found (for backwards compat)
            content.append(f"{idx} {coord[1]} {coord[2]} {coord[3]}")
        content.append("</pos>")
        return content
    
    
    def _add_scores_to_content(self, content):
        content.append("<hsp>")
        for index, row in self.scores.iterrows():
            # Convert UIDs back to indices for CLANS file format
            uid1 = row['PDBchain1']
            uid2 = row['PDBchain2']
            idx1 = self._uid_to_index.get(uid1, uid1)  # Fallback for backwards compat
            idx2 = self._uid_to_index.get(uid2, uid2)
            content.append(f"{idx1} {idx2}:{row['score']}")
        content.append("</hsp>")
        return content
