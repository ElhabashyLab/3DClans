import enum


class InputFileType(enum.Enum):
    """
    This enum class represents the possible types of a input file.
    """
    CLANS = "clans" 
    FASTA = "fasta" 
    A2M = "a2m" # The input file is in the format of an A2M file (treated the same as fasta)
    TSV = "tsv" # The input files is a tsv file containg the columns [uid, region_start, region_end, ...]
    CONF = "conf" # see ConfigFile.py class for more info
