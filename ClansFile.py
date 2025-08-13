class ClansFile:
    """
    This Class represents a CLANS file
    A Clans file consists of:
    - number of sequences
    - the fasta of the sequences
    - coordinates of the sequences (at first randomly generated)
    - pairwise similarity scores of the sequences
    
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
    def __init__(self, number_of_sequences, path_to_fasta, coordinates, scores):
        self.number_of_sequences = number_of_sequences
        self.path_to_fasta = path_to_fasta
        self.coordinates = coordinates
        self.scores = scores
        
        
    def __str__(self):
        """
        Returns a string representation of the CLANS file.
        """
        output = []
        # add the number of sequences
        output.append(f"sequences={self.number_of_sequences}")
        output.append("<seq>")
        # add the content of the fasta file
        with open(self.path_to_fasta, 'r') as fasta_file:
            for line in fasta_file:
                output.append(line.strip()) 
        output.append("</seq>")
        # add coordinates for each sequence
        output.append("<pos>")
        for coord in self.coordinates:
            output.append(f"{coord[0]} {coord[1]} {coord[2]} {coord[3]}")
        output.append("</pos>")
        # add pairwise similarity scores (TM) (scores is a df with columns = ['PDBchain1', 'PDBchain2', 'TM1', 'TM2', 'TM])
        
        # need to check if scores are in correct order
        
        output.append("<hsp>")
        for index, row in self.scores.iterrows():
            output.append(f"{row['PDBchain1']} {row['PDBchain2']}:{row['TM']}")     
        output.append("</hsp>")
        return "\n".join(output)
        