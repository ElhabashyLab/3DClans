import enum


class InputFileType(enum.Enum):
    """
    This enum class represents the possible types of a input file.
    """
    CLANS = "clans" # This is input type is used to run the recovered clans.jar in old_clans folder
    FASTA = "fasta" # The input file is a normal fasta file with all the sequences which should be compared
    A2M = "a2m" # The input file is in the format of an A2M files and contains already cut regions of the sequence
    TSV = "tsv" # The input files is a tsv file containg the columns [uid, region_start, region_end, ...]
    CONFIG = "config" # The input file is a config file containg all the necassary information
    