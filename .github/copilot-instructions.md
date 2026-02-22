# Clans-3D AI Coding Agent Instructions

## Project Overview

Clans-3D is a bioinformatics tool that generates CLANS visualization files from protein structures. It computes pairwise structural similarity scores and performs 2D/3D clustering for visualization with the CLANS software.

**Key Data Flow:** Input file → Protein extraction → Structure retrieval → Similarity computation → CLANS file generation

## Package Structure

The project uses a `src/clans3d/` package layout with `pyproject.toml` for packaging.

```
src/clans3d/
├── __init__.py
├── main.py                          # CLI entry point & create_clans_file()
├── core/                            # Core data types and file handling
│   ├── input_file_type.py           # InputFileType enum (FASTA, A2M, TSV, CLANS, CONF)
│   ├── config_file.py               # ConfigFile for -key value config format
│   ├── clans_file.py                # ClansFile data class
│   └── clans_file_generator.py      # Parse/generate CLANS files
├── similarity/                      # Structural similarity computation
│   ├── tool_type.py                 # ToolType enum (Foldseek, USalign)
│   ├── struct_sim_tool.py           # Abstract base class for tools
│   ├── foldseek.py                  # Foldseek implementation
│   ├── usalign.py                   # USalign implementation
│   └── struct_sim_computer.py       # Orchestrates tool selection/execution
├── utils/                           # Shared utilities
│   ├── file_utils.py                # File/directory operations
│   ├── api_utils.py                 # UniProt/UniParc API calls
│   ├── fasta_utils.py               # FASTA parsing and generation
│   ├── structure_utils.py           # PDB structure retrieval
│   └── dependency_checks.py         # External tool verification
├── evaluation/                      # Analysis and clustering
│   ├── scores_evaluator.py          # Main evaluation orchestrator
│   ├── clans_data_extractor.py      # Extract data from CLANS files
│   ├── clans_visualizer.py          # Scatter/bar/heatmap plots
│   ├── cluster_analyzer.py          # DBSCAN/HDBSCAN/Leiden clustering
│   └── data_normalizer.py           # min-max, z-score, -log10
├── legacy/                          # Legacy CLANS Java integration
│   ├── utils_old_clans.py           # Headless CLANS execution, BLAST-based CLANS
│   ├── clans_website_version.jar
│   └── clans_working_version.jar
├── dataset_generator/               # Synthetic test data generation
│   └── dataset_generator.py         # Generate clustered FASTA datasets via PSI-BLAST
└── benchmark/                       # Performance benchmarking
    └── benchmark_tool_speed.py      # Benchmark similarity tools
```

## Architecture & Core Components

### Input Pipeline

- **InputFileType enum** [`src/clans3d/core/input_file_type.py`]: Defines supported formats (FASTA, A2M, TSV, CLANS, CONF)
- **Supported input formats**:
  - **FASTA/A2M**: Protein sequences with region specification format `>id/start-end` (e.g., `>sp|P49811|MYOD1_PIG/2-300`)
  - **TSV**: Columns `[entry, region_start, region_end]`
  - **CLANS**: For loading existing CLANS files for re-analysis
  - **CONF**: Configuration files for batch processing

### Similarity Computation

- **ToolType enum** [`src/clans3d/similarity/tool_type.py`]: Foldseek, USalign
- **StructSimTool** [`src/clans3d/similarity/struct_sim_tool.py`]: Abstract base class with subclass pattern
  - Subclasses override `start_run()` for tool-specific initialization and `_parse_output()` for score extraction
  - All tools compute all pairwise scores (expected count: `n*(n-1)/2` where n = number of structures)
  - Returns DataFrame with columns `[seq_id_1, seq_id_2, score]`
- **StructSimComputer** [`src/clans3d/similarity/struct_sim_computer.py`]: Orchestrates tool selection and execution
  - Foldseek score parameter: `evalue` (default) or `TM` score
  - Validates score count and warns on mismatches

### CLANS File Format

- **ClansFile** [`src/clans3d/core/clans_file.py`]: Represents CLANS file structure
  - Sections: `<param>`, `<seq>`, `<pos>` (coordinates), `<hsp>` (pairwise scores)
  - Format: `seq_id1 seq_id2:score` (e.g., `0 1:4.6e-50`)
- **ClansFileGenerator** [`src/clans3d/core/clans_file_generator.py`]: Bidirectional parsing and generation
  - `parse_clans_file()`: Extracts all data into ClansFile object
  - `generate_clans_file_*()`: Creates CLANS files with coordinates and scores

### Evaluation & Analysis

- **ScoresEvaluator** [`src/clans3d/evaluation/scores_evaluator.py`]: Main orchestrator with component access
  - `extract_data_from_clans_files()`: Merges structural and sequence scores
  - Exposes `extractor`, `normalizer`, `clustering`, `visualizer` as public attributes
- **ClusterAnalyzer** [`src/clans3d/evaluation/cluster_analyzer.py`]: Clustering methods
  - `find_clusters_density_based()`: DBSCAN/HDBSCAN on 3D coordinates
  - `find_clusters_graph_based()`: Leiden community detection on score graphs
- **DataNormalizer** [`src/clans3d/evaluation/data_normalizer.py`]: min-max, z-score, -log10
- **ClansVisualizer** [`src/clans3d/evaluation/clans_visualizer.py`]: scatter, bar, heatmap plots

## Import Conventions

All imports use the `clans3d` package path:

```python
from clans3d.core.input_file_type import InputFileType
from clans3d.similarity.tool_type import ToolType
from clans3d.similarity.struct_sim_computer import StructSimComputer
from clans3d.utils.structure_utils import fetch_pdbs
from clans3d.evaluation.scores_evaluator import ScoresEvaluator
```

## Critical Workflows

### Main Workflow (main.py)

1. Parse CLI args or config file (ConfigFile handles both)
2. Call `create_clans_file()` which:
   - Validates and normalizes input file
   - Extracts protein UIDs/regions
   - Downloads/retrieves PDB structures
   - Runs similarity tool via StructSimComputer
   - Generates ClansFile with coordinates (from CLANS visualization)
3. Output: CLANS file + cleaned FASTA file

### CLI Entry Point

```bash
clans3d -l input.fasta -i fasta -t foldseek [-s evalue|TM] [-c config.conf]
```

### External Tool Integration

- **Subprocess execution** [`StructSimTool._execute_run()`]:
  - Tools run as subprocesses with captured output
  - Must handle tool-specific output parsing (Foldseek: TSV, USalign: custom format)
  - Tools expected in system PATH
- **Database/Index creation**:
  - Foldseek: Creates `.dbtype`, `.index`, `.lookup` files in working_dir
  - USalign: No pre-processing needed

### Configuration System

- **ConfigFile** [`src/clans3d/core/config_file.py`]: Reads/writes config with format `-key value`
- Config arguments merged with CLI args (CLI takes precedence)
- Used by CLANS headless execution and batch processing

## Key Patterns & Conventions

### Enum-Based Type Safety

- Use `ToolType` and `InputFileType` enums throughout, never string comparisons
- Convert from string via `ToolType(string_value)` and check `.value` property

### DataFrame-Centric Data Flow

- Pairwise scores always: `[Sequence_ID_1, Sequence_ID_2, Score]` columns
- Coordinates always: `[Sequence_ID, x, y, z]` columns
- Operations: merge on shared ID columns, normalize scores before analysis

### Suffix Conventions for Multi-Source Data

- `_struct`: Structural-based scores/coordinates
- `_seq`: Sequence-based scores (from BLAST in legacy)
- `_DBSCAN`, `_HDBSCAN`, `_Leiden`: Algorithm-specific clustering results

### Error Handling

- Warn on score count mismatches (don't fail; may indicate edge case)
- Validate input file format before processing
- Raise ValueError for missing required columns or configuration

### Working Directories

- Tool working directories cleaned/reset per run (see `reset_dir_content()` in utils)
- Coordinates computed from CLANS layout; not user-provided

### File Naming

- All Python source files use `snake_case` naming
- Module-specific enums live in their own files: `input_file_type.py` in core, `tool_type.py` in similarity

## Dependencies & External Requirements

- **Python 3.10+** with pip packages (see pyproject.toml): pandas, scipy, scikit-learn, hdbscan, networkx, igraph, leidenalg, biopython, matplotlib, seaborn
- **External tools** must be in PATH:
  - foldseek (fast similarity computation)
  - USalign (slower, structure alignment)
- **CLANS software**: For visualization (headless execution via `clans3d.legacy.utils_old_clans`)
- **PDB database access**: For structure retrieval (AlphaFold DB, online PDB queries)

## Testing Considerations

- Example files in `example_files/` (small_fasta_files, big_fasta_files, config_files)
- Test files in `clans_files/` (100.clans, 50_cleaned.clans, 500_cleaned.clans)
- Test scaffolding in `tests/` with `tests/fixtures/`
- Dataset generator: `clans3d.dataset_generator.dataset_generator.DatasetGenerator`
- Benchmark: `clans3d.benchmark.benchmark_tool_speed.Benchmark`

## Common Developer Tasks

1. **Adding a new similarity tool**: Extend `StructSimTool` in `src/clans3d/similarity/`, implement `start_run()` and `_parse_output()`, register in `StructSimComputer._set_up_tools()`
2. **Adding input format**: Add to `InputFileType` enum in `src/clans3d/core/input_file_type.py`, handle parsing in `src/clans3d/utils/fasta_utils.py` and `structure_utils.py`
3. **Modifying clustering**: Update `ClusterAnalyzer` methods in `src/clans3d/evaluation/cluster_analyzer.py`
4. **Debugging score discrepancies**: Check score count validation in `StructSimComputer` and DataFrame merge logic in `ScoresEvaluator.extract_data_from_clans_files()`
