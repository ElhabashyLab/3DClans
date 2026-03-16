# Clans-3D

## Table of Contents

- [How It Works](#how-it-works)
- [Installation](#installation)
  - [Prerequisites](#prerequisites)
  - [Installing Foldseek](#installing-foldseek)
  - [Installing USalign](#installing-usalign)
  - [Adding tools to your PATH](#adding-tools-to-your-path)
  - [Install from source](#install-from-source)
  - [Verify installation](#verify-installation)
- [Usage](#usage)
  - [Command-Line Interface](#command-line-interface)
  - [Input File Formats](#input-file-formats)
  - [Configuration Files](#configuration-files)
  - [Examples](#examples)
  - [Output](#output)
  - [CLANS File Format](#clans-file-format)
- [Benchmark](#benchmark)
- [Evaluation](#evaluation)
- [Legacy CLANS Integration](#legacy-clans-integration)
- [Project Structure](#project-structure)
- [Testing](#testing)
- [License](#license)

---

**Clans-3D** generates [CLANS](https://www.eb.tuebingen.mpg.de/protein-evolution/software/clans/) visualization files from protein structures using pairwise structural similarity scores. Instead of relying on sequence-based BLAST scores like the original CLANS tool, Clans-3D computes structural similarities with [Foldseek](https://github.com/steineggerlab/foldseek) or [USalign](https://zhanggroup.org/US-align/), enabling structure-aware clustering and visualization of protein families.

## How It Works

Clans-3D automates the full pipeline from a list of proteins to a ready-to-visualize CLANS file:

```
Input File (FASTA / A2M / A3M / TSV)
        │
        ▼
  Parse protein identifiers & regions
        │
        ▼
  Download AlphaFold structures (CIF format)
        │
        ▼
  Compute all-vs-all structural similarity
  (Foldseek or USalign)
        │
        ▼
  Generate CLANS file with scores & random 3D coordinates
        │
        ▼
  Output: .clans file  +  cleaned FASTA file
```

1. **Input parsing** — Extracts UniProt accessions and optional region annotations from the input file.
2. **Structure retrieval** — Downloads AlphaFold-predicted structures in CIF format. If a region is specified (e.g., residues 2–300), only that region is extracted. Downloads run in parallel (default: 10 workers), configurable via `-w/--workers`.
3. **Pairwise similarity** — Runs an all-vs-all structural comparison using Foldseek (fast, E-value or TM-score) or USalign (slower, TM-score). This produces `n*(n-1)/2` pairwise scores.
4. **CLANS file generation** — Writes a CLANS-format file containing the sequences, pairwise scores, and initial random 3D coordinates for visualization.

The output `.clans` file can be opened directly with one of the bundled CLANS JAR files in `src/clans3d/legacy/` (for example, `clans_working_version.jar`) for interactive 2D/3D clustering.

---

## Installation

### Prerequisites

- **Python ≥ 3.10**
- **Operating system**: Linux or macOS recommended. Windows is only partially supported — [USalign](https://zhanggroup.org/US-align/) provides a Windows binary, but [Foldseek](https://github.com/steineggerlab/foldseek) is **not available on Windows** (Linux and macOS only).
- At least one structural similarity tool in your system's `PATH`:
  - Foldseek
  - USalign

### Installing Foldseek

Download with:

| Platform      | Command                                                                                                    |
| ------------- | ---------------------------------------------------------------------------------------------------------- |
| Linux (AVX2)  | `wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz && tar xvzf foldseek-linux-avx2.tar.gz`       |
| Linux (ARM64) | `wget https://mmseqs.com/foldseek/foldseek-linux-arm64.tar.gz && tar xvzf foldseek-linux-arm64.tar.gz`     |
| macOS         | `wget https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz && tar xvzf foldseek-osx-universal.tar.gz` |

other foldseek binaries can be found at [https://mmseqs.com/foldseek](https://mmseqs.com/foldseek).

### Installing USalign

Download from [https://zhanggroup.org/US-align/](https://zhanggroup.org/US-align/).

**Compile from source (recommended, Linux & macOS):**

```bash
wget https://zhanggroup.org/US-align/bin/module/USalign.cpp
g++ -static -O3 -ffast-math -o USalign USalign.cpp   # Remove -static on macOS
```

**Or download a precompiled binary:**

| Platform | Download                                                                            |
| -------- | ----------------------------------------------------------------------------------- |
| Linux    | [USalignLinux64.zip](https://zhanggroup.org/US-align/bin/module/USalignLinux64.zip) |
| macOS    | [USalignMac.tar.gz](https://zhanggroup.org/US-align/bin/module/USalignMac.tar.gz)   |
| Windows  | [USalignWin64.zip](https://zhanggroup.org/US-align/bin/module/USalignWin64.zip)     |

### Adding tools to your PATH

Clans-3D needs `foldseek` and/or `USalign` to be on your system PATH:

**Linux / macOS:**

```bash
# Option 1: copy the binary to a standard location
sudo cp foldseek /usr/local/bin/
sudo cp USalign /usr/local/bin/

# Option 2: add the tool's directory to your shell profile permanently
# For bash, add to ~/.bashrc. For zsh, add to ~/.zshrc.
echo 'export PATH="/path/to/tool/directory:$PATH"' >> ~/.bashrc
source ~/.bashrc
```

**Windows (USalign only):**

1. Press `Win+R`, type `sysdm.cpl`, press Enter
2. Go to **Advanced → Environment Variables**
3. Under _User variables_, select **Path** and click **Edit**
4. Click **New** and add the folder containing `USalign.exe`
5. Click OK and restart your terminal

### Install Clans-3D from source

```bash
git clone https://github.com/ElhabashyLab/Clans-3D.git
cd Clans-3D

python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate

pip install .
```

This installs the `clans3d` command-line tool. To also install development/test dependencies:

```bash
pip install -e ".[dev]"
```

### Verify installation

```bash
# Check that the CLI is available
clans3d -h

# Check that external tools are found
which foldseek

which USalign # Linux/macOS
Get-Command USalign # Windows (PowerShell)
```

---

## Usage

### Command-Line Interface

```
clans3d -l <INPUT_FILE> -i <INPUT_TYPE> -t <TOOL> [OPTIONS]
```

#### Required arguments

| Argument                  | Description                                   |
| ------------------------- | --------------------------------------------- |
| `-l, --load <PATH>`       | Path to the input file                        |
| `-i, --input_type <TYPE>` | Input format: `fasta`, `a2m`, `a3m`, or `tsv` |
| `-t, --tool <TOOL>`       | Similarity tool: `foldseek` or `USalign`      |

#### Optional arguments

| Argument              | Description                                                                    |
| --------------------- | ------------------------------------------------------------------------------ |
| `-o, --out <PATH>`    | Output path for the CLANS file (file path or directory; default: `output/clans_files/`) |
| `-s, --score <SCORE>` | Foldseek score type: `evalue` (default) or `TM`. Only valid with `-t foldseek` |
| `-c, --conf <PATH>`   | Configuration file (CLI arguments override config values)                      |
| `-w, --workers <N>`   | Number of parallel threads for structure downloads (default: `10`)             |
| `-v, --verbose`       | Enable debug-level logging                                                     |
| `-q, --quiet`         | Suppress all output except errors                                              |

### Input File Formats

#### FASTA / A2M / A3M

Standard sequence files. The header can optionally specify a region of interest:

```
>sp|P49811|MYOD1_PIG/2-300 Myoblast determination protein 1
MELLSPPLRDVDLTGPDGSLCNFATADDFYDDPCF...
>sp|Q02346|MYOD1_RAT
MELLSPPLRDTDLLGPDGSLCSFATRDDFYDDPCF...
```

- Region format: `>identifier/start-end` (1-based, inclusive)
- If no region is given, the full-length structure is used
- A2M/A3M gap characters (`.` and `-`) are automatically stripped

#### TSV

Tab-separated file with protein accessions and optional regions:

```
entry       region_start    region_end
P49811      2               300
C5BN68      528             601
Q02346
```

- `region_start` and `region_end` are optional but must appear together if used
- Rows without regions use the full-length structure

### Configuration Files

Configuration files use the format `-key value`, one per line:

```conf
# Example config file
-load examples/small_fasta_files/5.fasta
-input_type fasta
-tool foldseek
-score evalue
-workers 20
-out results/my_output.clans
-verbose
```

Use with: `clans3d -c config.conf`

CLI arguments override any values set in the config file.

### Examples

#### Basic usage with Foldseek (E-value scores)

```bash
clans3d -l proteins.fasta -i fasta -t foldseek
```

#### Foldseek with TM-score

```bash
clans3d -l proteins.fasta -i fasta -t foldseek -s TM
```

#### USalign

```bash
clans3d -l proteins.fasta -i fasta -t USalign
```

#### A2M alignment input

```bash
clans3d -l alignment.a2m -i a2m -t foldseek
```

#### TSV input

```bash
clans3d -l proteins.tsv -i tsv -t foldseek -s evalue
```

#### Verbose output for debugging

```bash
clans3d -l proteins.fasta -i fasta -t foldseek -v
```

#### Increase parallel structure downloads

```bash
clans3d -l proteins.fasta -i fasta -t foldseek -w 30
```

Use higher values on fast networks for large inputs. If the AlphaFold API rate-limits requests, reduce the value.

#### Custom output path

```bash
# Write to a specific file
clans3d -l proteins.fasta -i fasta -t foldseek -o results/my_output.clans

# Write to a directory (filename auto-derived from input)
clans3d -l proteins.fasta -i fasta -t foldseek -o results/
```

### Output

After a successful run, the following files are produced:

```
output/clans_files/
└── <input_name>_cleaned.clans     # CLANS visualization file

work/
├── cleaned_input_storage/
│   └── <input_name>_cleaned.fasta # FASTA with only successfully downloaded structures
├── structures/                    # Downloaded AlphaFold CIF files
│   ├── P49811.cif
│   ├── Q02346.cif
│   └── ...
└── foldseek/  (or usalign/)       # Tool working files
```

Use `-o / --out` to change the output location. A path ending in `.clans` is treated as a file path; any other path is treated as a directory (the filename is derived from the input file name).

The cleaned FASTA file contains only the proteins for which structures were successfully downloaded. The CLANS file contains all pairwise structural similarity scores and can be opened in the CLANS Java application for interactive clustering.

### CLANS File Format

The generated `.clans` file follows the standard CLANS format:

```
sequences=4
<param>
</param>
<seq>
>sp|P29331|MYOD1_SHEEP/1-200 Myoblast determination protein 1
MELLSPPLRDVDLTGPDGSLCNFATADDFYDDPCF...
>sp|P49811|MYOD1_PIG/2-300 Myoblast determination protein 1
MELLSPPLRDVDLTGPDGSLCNFATADDFYDDPCF...
</seq>
<pos>
0 0.04 0.8 0.203
1 0.779 0.186 0.182
</pos>
<hsp>
0 1:1.487e-31
0 2:1.521e-29
1 2:1.531e-42
</hsp>
```

- `<param>` — CLANS parameters (initially empty)
- `<seq>` — FASTA sequences
- `<pos>` — 3D coordinates (initially random; refined by CLANS clustering)
- `<hsp>` — Pairwise similarity scores (`index1 index2:score`)

---

## Benchmark

The benchmark module measures the performance of each structural similarity tool (Foldseek with E-value, Foldseek with TM-score, USalign) across the full pipeline: structure download, score computation, and CLANS file generation.

**Quick start:**

```python
from clans3d.benchmark import Benchmark
from clans3d.core.input_file_type import InputFileType

benchmark = Benchmark("examples/small_fasta_files/50.fasta", InputFileType.FASTA)
benchmark.run_all_tools()
benchmark.print_results()
benchmark.export_csv()
```

Or from the command line:

```bash
python -m clans3d.benchmark.benchmark
```

**Example output:**

```
================================================================================
BENCHMARK RESULTS
================================================================================
Input: 50.fasta (50 structures)
Expected pairwise scores: 1225
================================================================================

         Tool Score Type  Structures  Download  Computation  Generation   Total  Scores  Success
     Foldseek     evalue          50      2.34        11.87        0.78   14.99    1225     True
     Foldseek         TM          50      2.31        18.42        0.79   21.52    1225     True
      USalign          -          50      2.35       127.84        0.79  130.98    1225     True

--------------------------------------------------------------------------------
SUMMARY:
  Fastest tool: Foldseek (evalue) - 14.99s total
  Slowest tool: USalign - 130.98s total
--------------------------------------------------------------------------------
```

Results are saved to `benchmark_output/run_YYYYMMDD_HHMMSS/` including a CSV export and the generated CLANS files.

For further details, see [docs/BENCHMARK_USAGE.md](docs/BENCHMARK_USAGE.md).

---

## Evaluation

The `evaluation/` directory contains tools for analyzing and comparing CLANS files — assessing clustering quality, computing agreement metrics, and generating visualizations. It includes components for data extraction, normalization (min-max, z-score, -log10), clustering (DBSCAN, HDBSCAN, Leiden), evaluation metrics (ARI, NMI, Jaccard), and plotting (scatter, bar, heatmap).

The evaluation workflow is designed for use in Jupyter notebooks (see `evaluation/notebooks/evaluation.ipynb`) and is intended for future work on systematically benchmarking clustering quality across different tools and datasets.

---

## Legacy CLANS Integration

Clans-3D includes utilities for running the original CLANS Java application in headless mode, useful for comparing structural similarity results with traditional sequence-based (BLAST) clustering.

```python
from clans3d.core.config_file import ConfigFile
from clans3d.legacy.utils_old_clans import run_clans_headless

config = ConfigFile("config.conf")
output_path = run_clans_headless(config, "path/to/clans.jar")
```

The headless mode clusters a CLANS file for a given number of rounds and writes the result with updated coordinates.

---

## Project Structure

```
src/clans3d/
├── main.py                          # CLI entry point
├── core/                            # Input parsing, pipeline, CLANS file I/O
│   ├── cli.py                       # Argument parsing
│   ├── pipeline.py                  # ClansPipeline orchestration
│   ├── input_file_type.py           # InputFileType enum
│   ├── config_file.py               # Config file reader/writer
│   ├── clans_file.py                # ClansFile data class
│   └── clans_file_generator.py      # Parse/generate CLANS files
├── similarity/                      # Structural similarity tools
│   ├── tool_type.py                 # ToolType enum (Foldseek, USalign)
│   ├── struct_sim_tool.py           # Abstract base class
│   ├── foldseek.py                  # Foldseek implementation
│   ├── usalign.py                   # USalign implementation
│   └── struct_sim_computer.py       # Tool selection & execution
├── utils/                           # Shared utilities
│   ├── file_utils.py                # File/directory operations
│   ├── api_utils.py                 # UniProt/UniParc API calls
│   ├── fasta_utils.py               # FASTA parsing & generation
│   ├── structure_utils.py           # AlphaFold structure download
│   ├── dependency_checks.py         # External tool verification
│   └── log.py                       # Logging setup
├── benchmark/                       # Performance benchmarking
│   ├── benchmark.py                 # Benchmark runner
│   └── benchmark_result.py          # Result data class
└── legacy/                          # CLANS Java headless integration
    ├── utils_old_clans.py           # Headless CLANS & BLAST utilities
    ├── CLANS_website_version.jar    # CLANS software from official website (buggy)
    └── CLANS.jar                    # CLANS software (not buggy)

evaluation/                          # Analysis & clustering (separate from main package)
├── evaluation_src/                  # Evaluation modules
└── notebooks/                       # Jupyter notebooks
```

---

## Testing

```bash
# Run all unit tests
pytest tests/unit/

# Run with coverage
pytest tests/unit/ --cov=src/clans3d --cov-report=term-missing

# Run integration tests (uses fixture files, no external tools needed)
pytest tests/integration/

# Run end-to-end tests (requires Foldseek/USalign in PATH)
pytest -m e2e
```

See [docs/TESTING_GUIDE.md](docs/TESTING_GUIDE.md) for full details on test organization and coverage.

---

## License

MIT
