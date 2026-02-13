#!/usr/bin/env python3
"""
Example script demonstrating how to use the Clans-3D benchmark.

This script can be run directly or modified to benchmark your own input files.
"""

from clans3d.benchmark import Benchmark
from clans3d.core.input_file_type import InputFileType


def run_benchmark_example():
    """Run a simple benchmark example with a small FASTA file."""
    
    # Example 1: Benchmark with FASTA file
    print("Example 1: Benchmarking with FASTA file")
    print("-" * 80)
    
    input_file = "examples/small_fasta_files/5.fasta"
    
    benchmark = Benchmark(input_file, InputFileType.FASTA)
    df = benchmark.run_all_tools()
    
    # Display results
    benchmark.print_results()
    
    # Export to CSV
    benchmark.export_csv()
    
    print("\n" + "=" * 80)
    print("Example 1 completed!")
    print("=" * 80 + "\n")
    
    # Example 2: Benchmark with TSV file (commented out - remove comments to run)
    # print("\nExample 2: Benchmarking with TSV file")
    # print("-" * 80)
    # 
    # input_file_tsv = "examples/small_tsv_files/100.tsv"
    # benchmark_tsv = Benchmark(input_file_tsv, InputFileType.TSV)
    # df_tsv = benchmark_tsv.run_all_tools()
    # benchmark_tsv.print_results()
    # benchmark_tsv.export_csv()


if __name__ == "__main__":
    run_benchmark_example()
