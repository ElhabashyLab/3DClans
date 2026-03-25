import enum


class TmMode(enum.Enum):
    """TM-score aggregation mode.

    The selected mode controls how ``TM1`` and ``TM2`` are combined before
    converting similarity to CLANS distance (``1 - TM``).
    """

    MIN = "min"
    MAX = "max"
    AVG = "mean"
