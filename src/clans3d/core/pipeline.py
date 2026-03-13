"""
ClansPipeline: Reusable pipeline for CLANS file generation.

Encapsulates the complete workflow:
1. Fetch protein structures
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
from clans3d.utils.structure_utils import fetch_structures
from clans3d.utils.fasta_utils import (
    copy_records_from_fasta,
    generate_fasta_from_uids_with_regions,
    generate_fasta_from_alignment_file,
)
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
        structures_dir: Directory for downloaded structures.
        output_dir: Directory for generated CLANS files.
        cleaned_input_storage: Directory for cleaned/intermediate input files.
        tool_working_dir: Base directory for similarity tool working files.
        download_workers: Maximum number of concurrent structure download threads. Default is 10.
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
        cleaned_input_storage: str = os.path.join("work", "cleaned_input_storage"),
        tool_working_dir: str = "work",
        download_workers: int = 10,
    ):
        self.input_file = input_file
        self.input_type = input_type
        self.tool = tool
        self.foldseek_score = foldseek_score
        self.verbose = verbose
        self.quiet = quiet
        self.structures_dir = structures_dir
        self.output_dir = output_dir
        self.cleaned_input_storage = cleaned_input_storage
        self.tool_working_dir = tool_working_dir
        self.download_workers = download_workers


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
        # validate input file exists
        if not os.path.exists(self.config.input_file):
            raise FileNotFoundError(f"Input file not found: {self.config.input_file}")
        
        # check if input_type is compatible with input file extension
        allowed_extensions = {
            InputFileType.FASTA: InputFileType.FASTA.value,
            InputFileType.A2M: InputFileType.A2M.value,
            InputFileType.A3M: InputFileType.A3M.value,
            InputFileType.TSV: InputFileType.TSV.value
            }
        ext = os.path.basename(self.config.input_file).split(".")[1].lower()  # get extension without dot
        if allowed_extensions[self.config.input_type] != ext:
            raise ValueError(
                f"Input file extension {ext} does not match expected "
                f"type {self.config.input_type.value}.")
        
        # check foldseek_score is only set when tool is foldseek
        if self.config.tool != ToolType.FOLDSEEK and self.config.foldseek_score is not None:
            raise ValueError(
                "foldseek_score can only be specified when tool is Foldseek."
            )

    # ------------------------------------------------------------------
    # Individual pipeline steps
    # ------------------------------------------------------------------

    def fetch_structures(self) -> dict[str, tuple[int, int] | None]:
        """Step 1: Download protein structures.

        Returns:
            dict mapping UID to optional (region_start, region_end).

        Raises:
            RuntimeError: If no structures could be downloaded.
        """
        logger.info("Fetching structures from AlphaFold...")
        os.makedirs(self.config.structures_dir, exist_ok=True)
        uids_with_regions = fetch_structures(
            self.config.input_file,
            self.config.input_type,
            self.config.structures_dir,
            max_workers=self.config.download_workers,
        )
        if not uids_with_regions:
            raise RuntimeError(
                "No structures could be downloaded. "
                "Check your input file and network access."
            )
        return uids_with_regions

    def generate_cleaned_fasta(
        self, uids_with_regions: dict[str, tuple[int, int] | None]
    ) -> str:
        """Step 2: Generate a cleaned FASTA file from the structures which were able to be downloaded.

        Args:
            uids_with_regions: Mapping uid to optional (region_start, region_end).

        Returns:
            Path to the cleaned FASTA file.
        """
        logger.info("Generating cleaned FASTA file with downloaded structures...")
        os.makedirs(self.config.cleaned_input_storage, exist_ok=True)
        input_file_name = os.path.basename(self.config.input_file).split(".")[0]
        cleaned_path = os.path.join(
            self.config.cleaned_input_storage, f"{input_file_name}_cleaned.fasta"
        )

        if self.config.input_type == InputFileType.FASTA:
            # Copy sequences directly from original FASTA
            uids = list(uids_with_regions.keys())
            return copy_records_from_fasta(self.config.input_file, uids, cleaned_path)
        elif self.config.input_type in (InputFileType.A2M, InputFileType.A3M):
            # Strip alignment gaps locally, keeping all residues
            return generate_fasta_from_alignment_file(
                uids_with_regions, cleaned_path, self.config.input_file
            )
        else:  # TSV - download sequences from UniProt
            return generate_fasta_from_uids_with_regions(
                uids_with_regions, cleaned_path
            )

    def compute_scores(self, structures) -> pd.DataFrame:
        """Step 3: Compute pairwise structural similarity scores.
        
        Args:
            structures: Directory containing the downloaded structure files.
        Returns:
            DataFrame with columns [Sequence_ID_1, Sequence_ID_2, score].
        """
        logger.info(f"Computing pairwise similarity scores with {self.config.tool.value}...")
        foldseek_score = self.config.foldseek_score or "evalue"
        computer = StructSimComputer(
            foldseek_score=foldseek_score,
            working_dir=self.config.tool_working_dir
        )
        scores = computer.run(self.config.tool, structures)
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
        scores = self.compute_scores(self.config.structures_dir)
        clans_file_path = self.generate_clans_file(scores, cleaned_fasta_path)
        return clans_file_path, cleaned_fasta_path
