"""
Benchmark module for Clans-3D structural similarity tools.

This module provides comprehensive benchmarking of the full Clans-3D pipeline,
measuring timing for PDB download, score computation, and CLANS file generation.

Supported tools:
- Foldseek (with evalue and TM score types)
- USalign

See `docs/BENCHMARK_USAGE.md` for detailed usage instructions and example results.
"""

import time
import os
from datetime import datetime
import pandas as pd

from clans3d.benchmark.benchmark_result import BenchmarkResult
from clans3d.core.pipeline import ClansPipeline, PipelineConfig
from clans3d.similarity.tool_type import ToolType
from clans3d.core.input_file_type import InputFileType


class Benchmark:
    """
    Benchmark class for measuring performance of Clans-3D structural similarity tools.
    
    This class wraps :class:`ClansPipeline` and times each pipeline step
    independently for detailed performance analysis.
    
    Structures are downloaded and the cleaned FASTA is generated once during
    initialization so that ``run_single_tool`` can be called directly without
    any manual preparation.
    """
    
    def __init__(self, input_file: str, input_file_type: InputFileType, 
                 output_dir: str = "benchmark_output"):
        """
        Initialize benchmark with input file.
        
        Downloads PDB structures and generates the cleaned FASTA immediately
        so that individual tool benchmarks can run without extra setup.
        
        Args:
            input_file: Path to input file (FASTA, A2M, or TSV format)
            input_file_type: Type of input file
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
        
        # Download structures and generate cleaned FASTA upfront
        self._init_pipeline = self._make_pipeline(ToolType.FOLDSEEK, None)
        
        start = time.perf_counter()
        self.uids_with_regions = self._init_pipeline.fetch_structures()
        self.num_structures = len([f for f in os.listdir(self.structures_dir) 
                                  if f.endswith('.pdb') or f.endswith('.cif')])
        self.cleaned_fasta_path = self._init_pipeline.generate_cleaned_fasta(
            self.uids_with_regions
        )
        self.pdb_download_time = time.perf_counter() - start
        print(f"Downloaded {self.num_structures} structures in {self.pdb_download_time:.2f}s")
    
    def _make_pipeline(self, tool: ToolType, score_type: str | None) -> ClansPipeline:
        """Create a ClansPipeline for the given tool configuration."""
        config = PipelineConfig(
            input_file=self.input_file,
            input_type=self.input_file_type,
            tool=tool,
            foldseek_score=score_type,
            structures_dir=self.structures_dir,
            output_dir=self.clans_dir,
            input_storage_dir=self.work_dir,
        )
        return ClansPipeline(config)
    
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
        
        # Run Foldseek with both score types
        print("Running similarity tools...")
        print("-" * 80)
        
        for score_type in ["evalue", "TM"]:
            result = self.run_single_tool(ToolType.FOLDSEEK, score_type)
            self.results.append(result)
        
        # Run USalign
        result = self.run_single_tool(ToolType.USALIGN, None)
        self.results.append(result)
        
        print("\n" + "="*80)
        print("Benchmark completed!")
        print("="*80 + "\n")
        
        return self.get_results_df()
    
    def run_single_tool(self, tool_type: ToolType, score_type: str | None) -> BenchmarkResult:
        """
        Run benchmark for a single tool configuration.
        
        Structures are already downloaded during initialization, so this
        only measures score computation and CLANS file generation.
        
        Args:
            tool_type: Type of tool to benchmark
            score_type: Score type for Foldseek ("evalue" or "TM"), None for other tools
            
        Returns:
            BenchmarkResult: Results from this tool run
        """
        tool_name = tool_type.value
        score_label = f" ({score_type})" if score_type else ""
        print(f"\nBenchmarking: {tool_name}{score_label}")
        
        pipeline = self._make_pipeline(tool_type, score_type)
        
        try:
            # Score Computation
            start = time.perf_counter()
            scores = pipeline.compute_scores()
            score_computation_time = time.perf_counter() - start
            print(f"Computed scores for {len(scores)} pairs in {score_computation_time:.2f}s")
            
            # CLANS Generation
            start = time.perf_counter()
            output_name = f"{tool_name}_{score_type if score_type else 'default'}.clans"
            pipeline.generate_clans_file(scores, self.cleaned_fasta_path,
                                         output_filename=output_name)
            clans_generation_time = time.perf_counter() - start
            print(f"Generated CLANS file {output_name} in {clans_generation_time:.2f}s")
            
            total_time = self.pdb_download_time + score_computation_time + clans_generation_time
            
            print(f"  - Total: {total_time:.2f}s | Scores: {len(scores)}")
            
            return BenchmarkResult(
                tool=tool_name,
                score_type=score_type,
                num_structures=self.num_structures,
                time_pdb_download=self.pdb_download_time,
                time_score_computation=score_computation_time,
                time_clans_generation=clans_generation_time,
                time_total=total_time,
                num_scores=len(scores),
                success=True,
                error_message=None
            )
            
        except Exception as e:
            print(f"Failed: {str(e)}")
            return BenchmarkResult(
                tool=tool_name,
                score_type=score_type,
                num_structures=self.num_structures,
                time_pdb_download=self.pdb_download_time,
                time_score_computation=0.0,
                time_clans_generation=0.0,
                time_total=self.pdb_download_time,
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
                "Success": "True" if result.success else "False"
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
    input_file = "examples/small_fasta_files/5.fasta"
    input_type = InputFileType.FASTA
    
    # Alternative TSV example:
    # input_file = "examples/small_tsv_files/5.tsv"
    # input_type = InputFileType.TSV
    
    # Run benchmark
    benchmark = Benchmark(input_file, input_type)
    df = benchmark.run_all_tools()
    
    # Display results
    benchmark.print_results()
    
    # Export to CSV
    benchmark.export_csv()
