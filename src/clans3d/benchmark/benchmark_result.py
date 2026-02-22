from dataclasses import dataclass

@dataclass
class BenchmarkResult:
    """Results from a single tool benchmark run."""
    tool: str                         # e.g., "Foldseek", "USalign"
    score_type: str | None            # "evalue" or "TM" for Foldseek, None otherwise
    num_structures: int               # Number of proteins processed
    time_structure_download: float     # Time to fetch structures (seconds)
    time_score_computation: float     # Time for similarity computation (seconds)
    time_clans_generation: float      # Time to generate CLANS file (seconds)
    time_total: float                 # Total end-to-end time (seconds)
    num_scores: int                   # Number of pairwise scores computed
    success: bool                     # Whether benchmark completed successfully
    error_message: str | None = None  # Error details if failed
    