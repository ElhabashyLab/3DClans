import enum


class InputFileType(enum.Enum):
    """
    This enum class represents the possible types of a input file.
    """
    CLANS = "clans" # This is input type is used to run clans.jar in the recovered_clans folder
    FASTA = "fasta" # The input file is a normal fasta file with all the sequences which should be compared
    A2M = "a2m" # The input file is in the format of an A2M file (treated the same as fasta)
    TSV = "tsv" # The input files is a tsv file containg the columns [uid, region_start, region_end, ...]
    CONF = "conf" # The input file is a config file containg all the necassary information (see ConfigFile.py class for more info)
    