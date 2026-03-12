import enum


class InputFileType(enum.Enum):
    """
    This enum class represents the possible types of a input file.
    """
    CLANS = "clans" 
    FASTA = "fasta" 
    A2M = "a2m"  # A2M alignment format
    A3M = "a3m"  # A3M alignment format (MMseqs2/ColabFold)
    TSV = "tsv"  # TSV file with columns [uid, region_start, region_end, ...]
    CONF = "conf"  # see ConfigFile.py class for more info
