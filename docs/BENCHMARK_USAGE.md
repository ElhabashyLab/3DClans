# Benchmark Usage Guide

## Overview

The Clans-3D benchmark provides comprehensive performance analysis of structural similarity tools across the full pipeline:

1. **Structure Download**: Time to fetch/retrieve protein structures
2. **Score Computation**: Time for pairwise similarity calculation
3. **CLANS Generation**: Time to create visualization file

Structures are downloaded automatically when the `Benchmark` is created, so individual tool runs require no manual preparation.

## Supported Tools

- **Foldseek** (with `evalue` and `TM` score types)
- **USalign**

## Usage

### Method 1: Use as Python Module

```python
from clans3d.benchmark import Benchmark
from clans3d.core.input_file_type import InputFileType

# Create benchmark instance (downloads structures automatically)
benchmark = Benchmark(
    input_file="path/to/input_file.fasta",
    input_file_type=InputFileType.FASTA
)

# Run all tools
df = benchmark.run_all_tools()

# Display results
benchmark.print_results()

# Export to CSV
benchmark.export_csv("my_benchmark_results.csv")
```

### Method 2: Run Directly

```bash
cd /path/to/Clans-3D
python -m clans3d.benchmark.benchmark
```

Note: This will perform `run_all_tools()`, `print_results()` and `export_csv()` with default output file name. Edit the `__main__` block in `src/clans3d/benchmark/benchmark.py` for changes.

## Input File Formats

Supported are FASTA, A2M, and TSV formats. See Input File Formats section in the main README for details.

## Output

### Console Output

The benchmark displays:

- Progress for each stage (structure download, computation, generation)
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
     Foldseek      evalue              50      2.34        11.87        0.78   14.99        1225       True
     Foldseek          TM              50      2.31        18.42        0.79   21.52        1225       True
      USalign           -              50      2.35       127.84        0.79  130.98        1225       True

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
    ├── structures/                          # Downloaded structure files
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

Structures are downloaded during `Benchmark.__init__()`, so `run_single_tool` works immediately:

```python
from clans3d.benchmark import Benchmark
from clans3d.core.input_file_type import InputFileType
from clans3d.similarity.tool_type import ToolType

# Initialize benchmark (downloads structures automatically)
benchmark = Benchmark("input.fasta", InputFileType.FASTA)

# Run single tool -- no manual download needed
result = benchmark.run_single_tool(ToolType.FOLDSEEK, "evalue")
print(f"Foldseek completed in {result.time_total:.2f}s")
```

## Additional Resources

- Main README: `/README.md`
- Benchmark source: `src/clans3d/benchmark/benchmark.py`
