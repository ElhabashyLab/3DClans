# Benchmark Restructure Plan (Simplified)

## Executive Summary

This document outlines a simplified plan to restructure the `benchmark_tool_speed.py` module into a professional benchmarking system for Clans-3D. The new benchmark will provide detailed timing diagnostics across all major pipeline stages with DataFrame-based output and optional CSV export.

---

## 1. Current State Analysis

### Existing Implementation Limitations

The current `benchmark_tool_speed.py` has several limitations:

1. **Limited timing granularity**: Only measures total tool execution time, not individual pipeline stages
2. **Confusing initialization**: The `run_with_structures_for_benchmark` flag creates unclear logic flow
3. **No CLANS file generation timing**: Doesn't measure CLANS file creation, which is part of the full workflow
4. **Incomplete workflow coverage**: Doesn't benchmark the full end-to-end pipeline (fetch → compute → generate)
5. **Poor result presentation**: Simple print statements, no structured output

### Current Workflow (from codebase analysis)

Based on `main.py::create_clans_file()`, the complete pipeline is:

```
Input File → fetch_pdbs() → StructSimComputer.run() → ClansFileGenerator.generate_clans_file() → Output
```

**Key timing points:**

- PDB structure download/retrieval (fetch_pdbs)
- Similarity score computation (StructSimComputer.run)
- CLANS file generation (ClansFileGenerator.generate_clans_file)

---

## 2. Design Goals (Simplified)

### Primary Objectives

1. **Comprehensive timing**: Measure all three major pipeline stages independently
2. **DataFrame output**: Results in pandas DataFrame for easy analysis
3. **Nice formatting**: Print results in readable table format
4. **CSV export**: Optional export to CSV file
5. **Flexible input**: Accept both FASTA and TSV formats
6. **All-tool comparison**: Run Foldseek (evalue/TM) and USalign automatically
7. **Simple interface**: Runnable Python script, no complex CLI

---

## 3. Proposed Architecture (Simplified)

### 3.1 Core Components

#### **BenchmarkResult** (Data Class)

Simple dataclass storing timing and metadata for a single benchmark run.

```python
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
```

#### **Benchmark** (Main Class)

Single class that handles everything: execution, timing, formatting, export.

**Key Methods:**

```python
class Benchmark:
    def __init__(self, input_file: str, input_file_type: InputFileType,
                 output_dir: str = "benchmark_output")

    def run_all_tools(self) -> pd.DataFrame
    def run_single_tool(self, tool_type: ToolType, score_type: str = None) -> BenchmarkResult
    def get_results_df(self) -> pd.DataFrame
    def print_results(self)
    def export_csv(self, output_path: str)
```

### 3.2 File Structure

Single file approach for simplicity:

```
src/clans3d/benchmark/
├── __init__.py
└── benchmark.py                 # BenchmarkResult dataclass + Benchmark class
```

---

## 4. Detailed Implementation Plan

### 4.1 Input Handling

**Supported formats:**

- FASTA files (protein sequences with optional region notation)
- TSV files (entry, region_start, region_end columns)

**Validation:**

- Check file exists and is readable
- Validate input type matches file format
- Ensure file contains at least 2 sequences (for pairwise comparison)

### 4.2 Timing Instrumentation Strategy

Use `time.perf_counter()` for high-resolution timing around each stage:

```python
import time

# Example timing wrapper
start = time.perf_counter()
result = function_to_time()
end = time.perf_counter()
elapsed = end - start
```

**Stage 1: PDB Download**

```python
def _time_pdb_fetch(self, input_file: str, input_type: InputFileType) -> tuple[dict, float]:
    """Returns (uids_with_regions, elapsed_time)"""
    start = time.perf_counter()
    uids_with_regions = fetch_pdbs(input_file, input_type, self.structures_dir)
    elapsed = time.perf_counter() - start
    return uids_with_regions, elapsed
```

**Stage 2: Score Computation**

```python
def _time_score_computation(self, tool_type: ToolType, pdb_dir: str,
                           score_type: str = None) -> tuple[DataFrame, float]:
    """Returns (scores_df, elapsed_time)"""
    computer = StructSimComputer(foldseek_score=score_type or "evalue")
    start = time.perf_counter()
    scores = computer.run(tool_type, pdb_dir)
    elapsed = time.perf_counter() - start
    return scores, elapsed
```

**Stage 3: CLANS Generation**

```python
def _time_clans_generation(self, scores: DataFrame, fasta_path: str,
                          out_path: str) -> tuple[str, float]:
    """Returns (clans_file_path, elapsed_time)"""
    generator = ClansFileGenerator()
    start = time.perf_counter()
    clans_path = generator.generate_clans_file(scores, fasta_path, out_path)
    elapsed = time.perf_counter() - start
    return clans_path, elapsed
```

### 4.3 Tool Iteration Strategy

**Tool Coverage:**

- Foldseek (with both `evalue` and `TM` score types) → 2 variations
- USalign → 1 variation
- **Total: 3 benchmark configurations**

**Note:** TMalign is excluded as it will be deprecated in future versions.

**Implementation:**

```python
def run_all_tools(self) -> pd.DataFrame:
    """Run benchmark for all tool configurations."""
    results = []

    # Foldseek with both score types
    for score_type in ["evalue", "TM"]:
        result = self.run_single_tool(ToolType.FOLDSEEK, score_type)
        results.append(result)

    # USalign
    result = self.run_single_tool(ToolType.USALIGN)
    results.append(result)

    return self.get_results_df()
```

### 4.4 Result Formatting

**DataFrame Output:**
Results stored in pandas DataFrame with columns:

- Tool
- Score Type
- Num Structures
- Download Time (s)
- Computation Time (s)
- Generation Time (s)
- Total Time (s)
- Num Scores
- Success

**Pretty Print Example:**

```
Benchmark Results
Input: 50.fasta (50 structures, 1225 expected scores)
================================================================================

         Tool Score Type  Num Structures  Download  Computation  Generation   Total  Num Scores  Success
     Foldseek      evalue              50      2.34        11.87        0.78   14.99        1225     True
     Foldseek          TM              50      2.31        18.42        0.79   21.52        1225     True
      USalign           -              50      2.35       127.84        0.79  130.98        1225     True

Summary:
- Fastest tool: Foldseek (evalue) - 14.99s total
- Slowest tool: USalign - 130.98s total
- Download time contributes ~2% to total time
```

### 4.5 Export Format

#### **CSV Export**

Simple CSV export using pandas `to_csv()`:

```csv
Tool,Score Type,Num Structures,Download Time (s),Computation Time (s),Generation Time (s),Total Time (s),Num Scores,Success
Foldseek,evalue,50,2.34,11.87,0.78,14.99,1225,True
Foldseek,TM,50,2.31,18.42,0.79,21.52,1225,True
USalign,-,50,2.35,127.84,0.79,130.98,1225,True
```

---

## 5. Usage Interface (Simplified)

### 5.1 As Python Script

Edit the `__main__` block in `benchmark.py`:

```python
if __name__ == "__main__":
    # Configuration
    input_file = "examples/small_fasta_files/50.fasta"
    input_type = InputFileType.FASTA

    # Run benchmark
    benchmark = Benchmark(input_file, input_type)
    df = benchmark.run_all_tools()

    # Display results
    benchmark.print_results()

    # Optional: Export to CSV
    benchmark.export_csv("benchmark_results.csv")
```

Run with:

```bash
python -m clans3d.benchmark.benchmark
# or
python src/clans3d/benchmark/benchmark.py
```

---

## 6. Integration with Existing Codebase

### 6.1 Dependencies on Existing Modules

The benchmark will reuse existing infrastructure:

**From `clans3d.utils.structure_utils`:**

- `fetch_pdbs()` - PDB download/retrieval

**From `clans3d.similarity`:**

- `StructSimComputer` - Similarity computation orchestration
- `ToolType` enum - Tool type definitions

**From `clans3d.core`:**

- `ClansFileGenerator` - CLANS file generation
- `InputFileType` enum - Input format definitions

**From `clans3d.utils.fasta_utils`:**

- `generate_fasta_from_uids_with_regions()` - Clean FASTA creation

**From `clans3d.utils.file_utils`:**

- `reset_dir_content()` - Directory cleanup

### 6.2 Minimal Code Changes Required

**No changes needed to existing modules!**

The benchmark will be a **pure consumer** of existing APIs, requiring zero modifications to core functionality. This is a key design principle for maintainability.

### 6.3 Working Directory Management

**Strategy:** Use isolated benchmark directories to avoid conflicts:

```
benchmark_output/
├── run_20260213_102345/
│   ├── structures/              # Downloaded PDB files
│   ├── work/                    # Tool working directories
│   ├── input_cleaned.fasta      # Cleaned FASTA
│   ├── results_foldseek_evalue.clans
│   ├── results_foldseek_TM.clans
│   ├── results_TMalign.clans
│   ├── results_USalign.clans
│   └── benchmark_report.json
└── run_20260213_145623/
    └── ...
```

This prevents contamination between runs and allows comparison of historical results.

---

## 7. Error Handling Strategy

### 7.1 Graceful Degradation

**Principle:** If one tool fails, continue benchmarking others.

```python
def run_single_tool(self, tool_type: ToolType, score_type: str = None) -> BenchmarkResult:
    try:
        # Execute benchmark stages
        ...
        return BenchmarkResult(success=True, ...)
    except Exception as e:
        return BenchmarkResult(
            success=False,
            error_message=str(e),
            time_total=0.0,
            ...
        )
```

---

## 8. Conclusion

This simplified plan provides a streamlined approach to creating a professional benchmarking tool. The design prioritizes:

- **Simplicity**: Single file, no CLI complexity
- **Comprehensive metrics**: Detailed timing breakdown
- **Flexibility**: Multiple input formats
- **Usability**: DataFrame output with nice printing and CSV export
- **Minimal invasiveness**: No changes to core codebase

The modular design allows for easy extension in the future if needed.

---

## 11. Implementation Checklist (Simplified)

### Phase 1: Core Implementation

- [ ] Create `BenchmarkResult` dataclass with timing fields
- [ ] Implement `Benchmark` class with timing wrappers
- [ ] Add support for Foldseek (evalue/TM) and USalign
- [ ] Add input validation (file exists, format matches)
- [ ] Implement error handling and graceful degradation
- [ ] Working directory management and cleanup

### Phase 2: Output & Export

- [ ] Implement `get_results_df()` to return pandas DataFrame
- [ ] Implement `print_results()` with nice formatting
- [ ] Add `export_csv()` functionality
- [ ] Test with example files

### Phase 3: Testing & Documentation

- [ ] Test with FASTA and TSV inputs
- [ ] Test error scenarios
- [ ] Update docstrings
- [ ] Add usage example in `__main__`

---

## 12. Success Criteria (Simplified)

The restructured benchmark will be considered successful when:

1. ✅ **Comprehensive timing**: All three stages measured independently
2. ✅ **Tools supported**: Foldseek (evalue/TM) and USalign (no TMalign)
3. ✅ **Input flexibility**: Accepts both FASTA and TSV formats
4. ✅ **DataFrame output**: Results as pandas DataFrame
5. ✅ **Nice formatting**: Clear printed table
6. ✅ **CSV export**: Optional CSV output
7. ✅ **Robust error handling**: Continues on individual tool failures
8. ✅ **Zero core changes**: Uses existing modules without modification
9. ✅ **Simple usage**: Easy to run and understand

---

## 13. Timeline Estimate (Simplified)

- **Core Implementation**: 2-3 hours
- **Output & Export**: 1 hour
- **Testing & Polish**: 1 hour
- **Total**: ~4-5 hours of focused development

---

## Conclusion

This plan provides a comprehensive roadmap for transforming the current benchmark into a professional, publication-ready tool. The design prioritizes:

- **Minimal invasiveness**: No changes to core codebase
- **Comprehensive metrics**: Detailed timing breakdown
- **Flexibility**: Multiple input formats and export options
- **Robustness**: Graceful error handling
- **Usability**: Simple CLI and clear output

The modular design allows for iterative implementation, starting with core functionality and progressively adding reporting and export features.
