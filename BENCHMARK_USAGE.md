# Benchmark Usage Guide

## Overview

The Clans-3D benchmark provides comprehensive performance analysis of structural similarity tools across the full pipeline:

1. **PDB Download**: Time to fetch/retrieve protein structures
2. **Score Computation**: Time for pairwise similarity calculation
3. **CLANS Generation**: Time to create visualization file

## Supported Tools

- **Foldseek** (with `evalue` and `TM` score types)
- **USalign**

Note: TMalign support has been removed as it will be deprecated in future versions.

## Quick Start

### Method 1: Run the Example Script

```bash
cd /path/to/Clans-3D
python example_benchmark.py
```

### Method 2: Use as Python Module

```python
from clans3d.benchmark import Benchmark
from clans3d.core.input_file_type import InputFileType

# Create benchmark instance
benchmark = Benchmark(
    input_file="examples/small_fasta_files/50.fasta",
    input_file_type=InputFileType.FASTA
)

# Run all tools
df = benchmark.run_all_tools()

# Display results
benchmark.print_results()

# Export to CSV
benchmark.export_csv("my_benchmark_results.csv")
```

### Method 3: Run Directly

```bash
cd /path/to/Clans-3D
python -m clans3d.benchmark.benchmark
```

Note: Edit the `__main__` block in `src/clans3d/benchmark/benchmark.py` to change the input file.

## Input File Formats

### FASTA Format

```fasta
>sp|P49811|MYOD1_PIG/2-300 Myogenic differentiation factor
MSPELLARLREELLRRAARHPEGKGRV...
>sp|P15172|MYOD1_MOUSE/2-318 Myogenic differentiation factor
MELLSPELLARLREELLRRAARHPEGK...
```

**Region specification:** Use format `>id/start-end` to specify protein regions.

### TSV Format

```tsv
entry	region_start	region_end
P49811	2	300
P15172	2	318
```

## Output

### Console Output

The benchmark displays:

- Progress for each stage (PDB download, computation, generation)
- Formatted results table with timing breakdown
- Summary statistics (fastest/slowest tool, timing breakdown)

Example:

```
================================================================================
BENCHMARK RESULTS
================================================================================
Input: 50.fasta (50 structures)
Expected pairwise scores: 1225
================================================================================

         Tool Score Type  Num Structures  Download  Computation  Generation   Total  Num Scores Success
     Foldseek      evalue              50      2.34        11.87        0.78   14.99        1225       ✓
     Foldseek          TM              50      2.31        18.42        0.79   21.52        1225       ✓
      USalign           -              50      2.35       127.84        0.79  130.98        1225       ✓

--------------------------------------------------------------------------------
SUMMARY:
  Fastest tool: Foldseek (evalue) - 14.99s total
  Slowest tool: USalign - 130.98s total
  Download time contributes ~2.1% to total time
--------------------------------------------------------------------------------
```

### CSV Export

Results are automatically saved as CSV in the benchmark output directory:

```
benchmark_output/run_20260213_150904/benchmark_results.csv
```

CSV columns:

- Tool
- Score Type
- Num Structures
- Download Time (s)
- Computation Time (s)
- Generation Time (s)
- Total Time (s)
- Num Scores
- Success

### Output Directory Structure

```
benchmark_output/
└── run_20260213_150904/
    ├── structures/                          # Downloaded PDB files
    ├── clans_files/                        # Generated CLANS files
    │   ├── foldseek_evalue.clans
    │   ├── foldseek_TM.clans
    │   └── USalign_default.clans
    ├── work/                               # Working files
    │   └── input_cleaned.fasta
    └── benchmark_results.csv               # CSV export
```

## Customization

### Change Input File

Edit the `__main__` block in `benchmark.py`:

```python
if __name__ == "__main__":
    # Change these values
    input_file = "path/to/your/input.fasta"
    input_type = InputFileType.FASTA  # or InputFileType.TSV

    benchmark = Benchmark(input_file, input_type)
    df = benchmark.run_all_tools()
    benchmark.print_results()
    benchmark.export_csv()
```

### Change Output Directory

```python
benchmark = Benchmark(
    input_file="input.fasta",
    input_file_type=InputFileType.FASTA,
    output_dir="my_custom_benchmark_dir"  # Custom output directory
)
```

### Run Single Tool

```python
from clans3d.similarity.tool_type import ToolType

# Initialize benchmark
benchmark = Benchmark("input.fasta", InputFileType.FASTA)

# Download PDBs first
benchmark.uids_with_regions = fetch_pdbs(
    benchmark.input_file,
    benchmark.input_file_type,
    benchmark.structures_dir
)

# Run single tool
result = benchmark.run_single_tool(ToolType.FOLDSEEK, "evalue", 0.0)
print(f"Foldseek completed in {result.time_total:.2f}s")
```

## Interpreting Results

### Timing Breakdown

- **Download Time**: Usually consistent across tools (same structures downloaded)
- **Computation Time**: Varies significantly:
  - Foldseek (evalue): Fastest, optimized for speed
  - Foldseek (TM): Slower than evalue, more detailed alignment
  - USalign: Slowest, most comprehensive structural alignment
- **Generation Time**: Usually negligible (<1s), just file I/O

### Performance Comparison

**For small datasets (<100 structures):**

- Foldseek (evalue): ~10-20x faster than USalign
- Foldseek (TM): ~5-10x faster than USalign

**For large datasets (>500 structures):**

- Speed difference becomes more pronounced
- Foldseek can complete in minutes vs hours for USalign

### When to Use Each Tool

- **Foldseek (evalue)**: Best for fast exploration, large datasets
- **Foldseek (TM)**: Balance of speed and alignment quality
- **USalign**: Most accurate, best for detailed analysis and publication

## Troubleshooting

### Tool Not Found

If you see errors like "foldseek: command not found":

1. Install the tool (see main README)
2. Ensure it's in your system PATH
3. Verify installation: `which foldseek`

### PDB Download Failures

- Check internet connection
- Some UniProt entries may not have structures
- The benchmark continues with successfully downloaded structures

### Out of Memory

For very large datasets (>1000 structures):

- Use a machine with more RAM
- Run tools individually instead of all at once
- Consider batching the input file

## Example Benchmarks

### Small Dataset (5 structures)

- Download: ~1-2s
- Foldseek: ~1-2s (computation)
- USalign: ~5-10s (computation)

### Medium Dataset (50 structures)

- Download: ~5-10s
- Foldseek: ~10-20s (computation)
- USalign: ~2-5 minutes (computation)

### Large Dataset (500 structures)

- Download: ~30-60s
- Foldseek: ~2-5 minutes (computation)
- USalign: ~1-3 hours (computation)

## Additional Resources

- Main README: `/README.md`
- Benchmark source: `src/clans3d/benchmark/benchmark.py`
- Example script: `example_benchmark.py`
- Architecture plan: `BENCHMARK_RESTRUCTURE_PLAN.md`
