"""
Benchmark module for Clans-3D structural similarity tools.

This module provides comprehensive benchmarking of the full Clans-3D pipeline,
measuring timing for PDB download, score computation, and CLANS file generation.

Supported tools:
- Foldseek (with evalue and TM score types)
- USalign

Usage:
    Edit the __main__ block with your input file and run:
    python -m clans3d.benchmark.benchmark
"""

import time
import os
from dataclasses import dataclass, asdict
from datetime import datetime
import pandas as pd

from clans3d.utils.structure_utils import fetch_pdbs
from clans3d.utils.fasta_utils import generate_fasta_from_uids_with_regions
from clans3d.similarity.struct_sim_computer import StructSimComputer
from clans3d.similarity.tool_type import ToolType
from clans3d.core.clans_file_generator import ClansFileGenerator
from clans3d.core.input_file_type import InputFileType


@dataclass
class BenchmarkResult:
    """Results from a single tool benchmark run."""
    tool: str                         # e.g., "Foldseek", "USalign"
    score_type: str | None            # "evalue" or "TM" for Foldseek, None otherwise
    num_structures: int               # Number of proteins processed
    time_pdb_download: float          # Time to fetch PDB structures (seconds)
    time_score_computation: float     # Time for similarity computation (seconds)
    time_clans_generation: float      # Time to generate CLANS file (seconds)
    time_total: float                 # Total end-to-end time (seconds)
    num_scores: int                   # Number of pairwise scores computed
    success: bool                     # Whether benchmark completed successfully
    error_message: str | None = None  # Error details if failed


class Benchmark:
    """
    Benchmark class for measuring performance of Clans-3D structural similarity tools.
    
    This class orchestrates the full pipeline:
    1. PDB structure download/retrieval
    2. Similarity score computation
    3. CLANS file generation
    
    Each stage is timed independently for detailed performance analysis.
    """
    
    def __init__(self, input_file: str, input_file_type: InputFileType, 
                 output_dir: str = "benchmark_output"):
        """
        Initialize benchmark with input file.
        
        Args:
            input_file: Path to input file (FASTA or TSV format)
            input_file_type: Type of input file (InputFileType.FASTA or InputFileType.TSV)
            output_dir: Directory for benchmark outputs (default: "benchmark_output")
        """
        if not os.path.exists(input_file):
            raise FileNotFoundError(f"Input file not found: {input_file}")
        
        self.input_file = input_file
        self.input_file_type = input_file_type
        self.output_dir = output_dir
        
        # Create timestamped output directory
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        self.run_dir = os.path.join(output_dir, f"run_{timestamp}")
        os.makedirs(self.run_dir, exist_ok=True)
        
        # Subdirectories for this benchmark run
        self.structures_dir = os.path.join(self.run_dir, "structures")
        self.clans_dir = os.path.join(self.run_dir, "clans_files")
        self.work_dir = os.path.join(self.run_dir, "work")
        
        os.makedirs(self.structures_dir, exist_ok=True)
        os.makedirs(self.clans_dir, exist_ok=True)
        os.makedirs(self.work_dir, exist_ok=True)
        
        # Store results
        self.results: list[BenchmarkResult] = []
        self.uids_with_regions = None
        self.cleaned_fasta_path = None
        self.num_structures = 0
    
    def run_all_tools(self) -> pd.DataFrame:
        """
        Run benchmark for all supported tool configurations.
        
        Returns:
            pd.DataFrame: Results as a pandas DataFrame
        """
        print(f"\n{'='*80}")
        print(f"Starting Clans-3D Benchmark")
        print(f"{'='*80}")
        print(f"Input file: {self.input_file}")
        print(f"Input type: {self.input_file_type.value}")
        print(f"Output directory: {self.run_dir}")
        print(f"{'='*80}\n")
        
        # Download PDBs once (shared across all tools)
        print("Stage 1/3: Downloading PDB structures...")
        start = time.perf_counter()
        try:
            self.uids_with_regions = fetch_pdbs(
                self.input_file, 
                self.input_file_type, 
                self.structures_dir
            )
            self.num_structures = len([f for f in os.listdir(self.structures_dir) 
                                      if f.endswith('.pdb') or f.endswith('.cif')])
            pdb_download_time = time.perf_counter() - start
            print(f"✓ Downloaded {self.num_structures} structures in {pdb_download_time:.2f}s\n")
        except Exception as e:
            print(f"✗ PDB download failed: {e}")
            return self.get_results_df()
        
        # Create cleaned FASTA file (needed for CLANS generation)
        input_file_name = os.path.basename(self.input_file).split(".")[0]
        self.cleaned_fasta_path = os.path.join(self.work_dir, f"{input_file_name}_cleaned.fasta")
        
        if self.input_file_type == InputFileType.FASTA:
            self.cleaned_fasta_path = generate_fasta_from_uids_with_regions(
                self.uids_with_regions, 
                self.cleaned_fasta_path, 
                self.input_file
            )
        else:
            self.cleaned_fasta_path = generate_fasta_from_uids_with_regions(
                self.uids_with_regions, 
                self.cleaned_fasta_path
            )
        
        # Run Foldseek with both score types
        print("Stage 2/3: Running similarity tools...")
        print("-" * 80)
        
        for score_type in ["evalue", "TM"]:
            result = self.run_single_tool(ToolType.FOLDSEEK, score_type, pdb_download_time)
            self.results.append(result)
        
        # Run USalign
        result = self.run_single_tool(ToolType.USALIGN, None, pdb_download_time)
        self.results.append(result)
        
        print("\n" + "="*80)
        print("Benchmark completed!")
        print("="*80 + "\n")
        
        return self.get_results_df()
    
    def run_single_tool(self, tool_type: ToolType, score_type: str | None, 
                       pdb_download_time: float) -> BenchmarkResult:
        """
        Run benchmark for a single tool configuration.
        
        Args:
            tool_type: Type of tool to benchmark
            score_type: Score type for Foldseek ("evalue" or "TM"), None for other tools
            pdb_download_time: Time already spent downloading PDBs (reused across tools)
            
        Returns:
            BenchmarkResult: Results from this tool run
        """
        tool_name = tool_type.value
        score_label = f" ({score_type})" if score_type else ""
        print(f"\nBenchmarking: {tool_name}{score_label}")
        
        try:
            # Stage 2: Score Computation
            print(f"  - Computing similarity scores...", end=" ", flush=True)
            start = time.perf_counter()
            
            computer = StructSimComputer(foldseek_score=score_type or "evalue")
            scores = computer.run(tool_type, self.structures_dir)
            
            score_computation_time = time.perf_counter() - start
            print(f"✓ ({score_computation_time:.2f}s)")
            
            if scores is None or len(scores) == 0:
                raise ValueError("No scores computed")
            
            # Stage 3: CLANS Generation
            print(f"  - Generating CLANS file...", end=" ", flush=True)
            start = time.perf_counter()
            
            generator = ClansFileGenerator()
            output_name = f"{tool_name}_{score_type if score_type else 'default'}.clans"
            out_path = os.path.join(self.clans_dir, output_name)
            
            if self.cleaned_fasta_path is None:
                raise ValueError("Cleaned FASTA path not initialized")
            
            clans_path = generator.generate_clans_file(scores, self.cleaned_fasta_path, out_path)
            
            clans_generation_time = time.perf_counter() - start
            print(f"✓ ({clans_generation_time:.2f}s)")
            
            total_time = pdb_download_time + score_computation_time + clans_generation_time
            
            print(f"  - Total: {total_time:.2f}s | Scores: {len(scores)}")
            
            return BenchmarkResult(
                tool=tool_name,
                score_type=score_type,
                num_structures=self.num_structures,
                time_pdb_download=pdb_download_time,
                time_score_computation=score_computation_time,
                time_clans_generation=clans_generation_time,
                time_total=total_time,
                num_scores=len(scores),
                success=True,
                error_message=None
            )
            
        except Exception as e:
            print(f"✗ Failed: {str(e)}")
            return BenchmarkResult(
                tool=tool_name,
                score_type=score_type,
                num_structures=self.num_structures,
                time_pdb_download=pdb_download_time,
                time_score_computation=0.0,
                time_clans_generation=0.0,
                time_total=pdb_download_time,
                num_scores=0,
                success=False,
                error_message=str(e)
            )
    
    def get_results_df(self) -> pd.DataFrame:
        """
        Convert results to pandas DataFrame.
        
        Returns:
            pd.DataFrame: Results with columns for each timing metric
        """
        if not self.results:
            return pd.DataFrame()
        
        data = []
        for result in self.results:
            row = {
                "Tool": result.tool,
                "Score Type": result.score_type if result.score_type else "-",
                "Num Structures": result.num_structures,
                "Download Time (s)": f"{result.time_pdb_download:.2f}",
                "Computation Time (s)": f"{result.time_score_computation:.2f}",
                "Generation Time (s)": f"{result.time_clans_generation:.2f}",
                "Total Time (s)": f"{result.time_total:.2f}",
                "Num Scores": result.num_scores,
                "Success": "✓" if result.success else "✗"
            }
            data.append(row)
        
        return pd.DataFrame(data)
    
    def print_results(self):
        """Print results in a nicely formatted table."""
        df = self.get_results_df()
        
        if df.empty:
            print("No results to display.")
            return
        
        print("\n" + "="*80)
        print("BENCHMARK RESULTS")
        print("="*80)
        print(f"Input: {os.path.basename(self.input_file)} ({self.num_structures} structures)")
        
        expected_scores = (self.num_structures * (self.num_structures - 1)) // 2
        print(f"Expected pairwise scores: {expected_scores}")
        print("="*80)
        print()
        
        # Print DataFrame with nice formatting
        pd.set_option('display.max_columns', None)
        pd.set_option('display.width', None)
        pd.set_option('display.max_colwidth', None)
        print(df.to_string(index=False))
        print()
        
        # Summary statistics
        successful_results = [r for r in self.results if r.success]
        if successful_results:
            fastest = min(successful_results, key=lambda r: r.time_total)
            slowest = max(successful_results, key=lambda r: r.time_total)
            
            print("-" * 80)
            print("SUMMARY:")
            print(f"  Fastest tool: {fastest.tool}" + 
                  (f" ({fastest.score_type})" if fastest.score_type else "") + 
                  f" - {fastest.time_total:.2f}s total")
            print(f"  Slowest tool: {slowest.tool}" + 
                  (f" ({slowest.score_type})" if slowest.score_type else "") + 
                  f" - {slowest.time_total:.2f}s total")
            
            avg_download = sum(r.time_pdb_download for r in successful_results) / len(successful_results)
            avg_total = sum(r.time_total for r in successful_results) / len(successful_results)
            download_percentage = (avg_download / avg_total) * 100 if avg_total > 0 else 0
            print(f"  Download time contributes ~{download_percentage:.1f}% to total time")
            print("-" * 80)
        
        print(f"\nResults saved to: {self.run_dir}")
        print("="*80 + "\n")
    
    def export_csv(self, output_path: str | None = None):
        """
        Export results to CSV file.
        
        Args:
            output_path: Path for CSV output. If None, saves to benchmark run directory.
        """
        df = self.get_results_df()
        
        if df.empty:
            print("No results to export.")
            return
        
        if output_path is None:
            output_path = os.path.join(self.run_dir, "benchmark_results.csv")
        
        df.to_csv(output_path, index=False)
        print(f"Results exported to: {output_path}")


if __name__ == "__main__":
    # Configuration - edit these values to benchmark different files
    input_file = "examples/small_fasta_files/50.fasta"
    input_type = InputFileType.FASTA
    
    # Alternative TSV example:
    # input_file = "examples/small_tsv_files/100.tsv"
    # input_type = InputFileType.TSV
    
    # Run benchmark
    benchmark = Benchmark(input_file, input_type)
    df = benchmark.run_all_tools()
    
    # Display results
    benchmark.print_results()
    
    # Export to CSV
    benchmark.export_csv()
