import enum


class InputFileType(enum.Enum):
    """
    This enum class represents the possible types of a input file.
    """
    CLANS = "clans" # This input type does not have any use yet
    FASTA = "fasta" # The input file is a norma fasta file with all the sequences which should be compared
    TSV = "tsv" # The input files is a tsv file containg the columns [uid, region_start, region_end, ...]
    CONFIG = "config" # The input file is a config file containg all the necassary information
    