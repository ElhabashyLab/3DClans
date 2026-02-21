"""
ClansPipeline: Reusable pipeline for CLANS file generation.

Encapsulates the complete workflow:
1. Fetch PDB structures
2. Generate cleaned FASTA
3. Compute pairwise similarity scores
4. Generate CLANS visualization file

Each step is exposed as a public method so callers (e.g. Benchmark)
can invoke and time them independently.
"""

import logging
import os

import pandas as pd

from clans3d.core.input_file_type import InputFileType
from clans3d.core.clans_file_generator import ClansFileGenerator
from clans3d.similarity.struct_sim_computer import StructSimComputer
from clans3d.similarity.tool_type import ToolType
from clans3d.utils.structure_utils import fetch_pdbs
from clans3d.utils.fasta_utils import generate_fasta_from_uids_with_regions
from clans3d.utils.log import setup_logging

logger = logging.getLogger(__name__)


class PipelineConfig:
    """Configuration for the CLANS generation pipeline.
    
    Args:
        input_file: Path to the input file (FASTA, A2M, or TSV).
        input_type: Format of the input file.
        tool: Structural similarity tool to use.
        foldseek_score: Score type for Foldseek ("evalue" or "TM"). Ignored for other tools.
        verbose: If ``True``, enable debug-level logging output.
        quiet: If ``True``, suppress all output except errors.  Overrides *verbose*.
        structures_dir: Directory for downloaded PDB structures.
        output_dir: Directory for generated CLANS files.
        input_storage_dir: Directory for cleaned/intermediate input files.
    """

    def __init__(
        self,
        input_file: str,
        input_type: InputFileType,
        tool: ToolType,
        foldseek_score: str | None = None,
        verbose: bool = False,
        quiet: bool = False,
        structures_dir: str = os.path.join("work", "structures"),
        output_dir: str = os.path.join("output", "clans_files"),
        input_storage_dir: str = os.path.join("work", "input_file_storage"),
    ):
        self.input_file = input_file
        self.input_type = input_type
        self.tool = tool
        self.foldseek_score = foldseek_score
        self.verbose = verbose
        self.quiet = quiet
        self.structures_dir = structures_dir
        self.output_dir = output_dir
        self.input_storage_dir = input_storage_dir


class ClansPipeline:
    """Orchestrates the CLANS file generation workflow.

    Usage as a library::

        config = PipelineConfig(
            input_file="proteins.fasta",
            input_type=InputFileType.FASTA,
            tool=ToolType.FOLDSEEK,
        )
        pipeline = ClansPipeline(config)
        clans_path, fasta_path = pipeline.run()

    Individual steps can also be called separately for benchmarking or
    custom workflows.
    """

    def __init__(self, config: PipelineConfig):
        self.config = config
        # Ensure logging is configured when used as a library (f.e. for benchmark)
        setup_logging(verbose=config.verbose, quiet=config.quiet)
        self._validate()

    def _validate(self) -> None:
        """Validate pipeline configuration."""
        if not os.path.exists(self.config.input_file):
            raise FileNotFoundError(
                f"Input file not found: {self.config.input_file}"
            )
        if (
            self.config.tool != ToolType.FOLDSEEK
            and self.config.foldseek_score is not None
        ):
            raise ValueError(
                "foldseek_score can only be specified when tool is Foldseek."
            )

    # ------------------------------------------------------------------
    # Individual pipeline steps
    # ------------------------------------------------------------------

    def fetch_structures(self) -> dict[str, tuple[int, int] | None]:
        """Step 1: Download PDB structures.

        Returns:
            dict mapping UID to optional (region_start, region_end).
        """
        logger.info("Fetching structures from AlphaFold...")
        os.makedirs(self.config.structures_dir, exist_ok=True)
        return fetch_pdbs(
            self.config.input_file,
            self.config.input_type,
            self.config.structures_dir,
        )

    def generate_cleaned_fasta(
        self, uids_with_regions: dict[str, tuple[int, int] | None]
    ) -> str:
        """Step 2: Generate a cleaned FASTA file from the fetched UIDs.

        Args:
            uids_with_regions: Mapping uid to optional (region_start, region_end).

        Returns:
            Path to the cleaned FASTA file.
        """
        os.makedirs(self.config.input_storage_dir, exist_ok=True)
        input_file_name = os.path.basename(self.config.input_file).split(".")[0]
        cleaned_path = os.path.join(
            self.config.input_storage_dir, f"{input_file_name}_cleaned.fasta"
        )

        if self.config.input_type == InputFileType.FASTA:
            return generate_fasta_from_uids_with_regions(
                uids_with_regions, cleaned_path, self.config.input_file
            )
        else:
            return generate_fasta_from_uids_with_regions(
                uids_with_regions, cleaned_path
            )

    def compute_scores(self) -> pd.DataFrame:
        """Step 3: Compute pairwise structural similarity scores.

        Returns:
            DataFrame with columns [PDBchain1, PDBchain2, score].
        """
        logger.info("Computing pairwise similarity scores with %s...", self.config.tool.value)
        foldseek_score = self.config.foldseek_score or "evalue"
        computer = StructSimComputer(foldseek_score=foldseek_score)
        scores = computer.run(self.config.tool, self.config.structures_dir)
        return scores

    def generate_clans_file(
        self, scores: pd.DataFrame, cleaned_fasta_path: str,
        output_filename: str | None = None,
    ) -> str:
        """Step 4: Generate a CLANS visualization file.

        Args:
            scores: Pairwise similarity scores.
            cleaned_fasta_path: Path to the cleaned FASTA.
            output_filename: Optional custom filename for the CLANS file.
                If ``None``, the name is derived from the cleaned FASTA basename.

        Returns:
            Path to the generated CLANS file.
        """
        logger.info("Generating CLANS file...")
        os.makedirs(self.config.output_dir, exist_ok=True)
        if output_filename is None:
            output_filename = (
                os.path.basename(cleaned_fasta_path).split(".")[0] + ".clans"
            )
        out_path = os.path.join(self.config.output_dir, output_filename)
        generator = ClansFileGenerator()
        path_to_clans_file = generator.generate_clans_file(scores, cleaned_fasta_path, out_path)
        logger.info("CLANS file generated at %s", path_to_clans_file)
        return path_to_clans_file

    # ------------------------------------------------------------------
    # Full pipeline execution
    # ------------------------------------------------------------------
    
    def run(self) -> tuple[str, str]:
        """Execute the full pipeline end-to-end.

        Returns:
            Tuple of (clans_file_path, cleaned_fasta_path).
        """
        uids_with_regions = self.fetch_structures()
        cleaned_fasta_path = self.generate_cleaned_fasta(uids_with_regions)
        scores = self.compute_scores()
        clans_file_path = self.generate_clans_file(scores, cleaned_fasta_path)
        return clans_file_path, cleaned_fasta_path
