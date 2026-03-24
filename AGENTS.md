# AGENTS.md — 3DClans Agent Guide

This file is the primary reference for AI coding agents working in this repository.
Read it **before** making any changes.

---

## Table of Contents

1. [What is 3DClans?](#1-what-is-3dclans)
2. [Repository Layout](#2-repository-layout)
3. [Environment Setup](#3-environment-setup)
4. [Running Tests](#4-running-tests)
5. [CLI Usage](#5-cli-usage)
6. [Architecture & Key Patterns](#6-architecture--key-patterns)
7. [Data Contracts](#7-data-contracts)
8. [Adding New Features](#8-adding-new-features)
9. [Code Conventions](#9-code-conventions)
10. [External Dependencies](#10-external-dependencies)
11. [Common Gotchas](#11-common-gotchas)

---

## 1. What is 3DClans?

3DClans generates [CLANS](https://www.eb.tuebingen.mpg.de/protein-evolution/software/clans/) visualization files from protein structures using **pairwise structural similarity** scores rather than sequence-based BLAST scores.

**Data flow:**

```
Input file (FASTA / A2M / A3M / TSV)
        │
        ▼
  Parse UniProt accessions & optional region annotations
        │
        ▼
  Download AlphaFold CIF structures (parallel, 10 workers default)
        │
        ▼
  All-vs-all structural similarity  ← Foldseek (evalue or TM) | USalign
        │
        ▼
  Generate .clans file  +  cleaned FASTA
```

The output `.clans` file can be opened in the bundled CLANS Java application for interactive 2D/3D clustering.

---

## 2. Repository Layout

```
3DClans/
├── AGENTS.md                         # ← you are here
├── README.md
├── pyproject.toml                    # packaging, deps, pytest config, coverage config
├── requirements.txt
├── docs/
│   ├── TESTING_GUIDE.md              # full test inventory with table per module
│   └── BENCHMARK_USAGE.md
├── examples/
│   ├── small_fasta_files/            # 5.fasta, 10.fasta, …  (real UniProt accessions)
│   ├── big_fasta_files/
│   ├── small_tsv_files/
│   ├── big_tsv_files/
│   ├── clans_files/                  # pre-generated .clans examples
│   ├── config_files/                 # example .conf files
│   └── generated_fasta/
├── src/clans3d/                      # installable package (pip install -e .)
│   ├── __init__.py
│   ├── main.py                       # CLI entry point → main()
│   ├── core/
│   │   ├── cli.py                    # argparse setup; parse_args()
│   │   ├── pipeline.py               # ClansPipeline + PipelineConfig
│   │   ├── input_file_type.py        # InputFileType enum (FASTA, A2M, A3M, TSV, CLANS, CONF)
│   │   ├── config_file.py            # ConfigFile: read/write -key value format
│   │   ├── clans_file.py             # ClansFile data class + __str__ serializer
│   │   └── clans_file_generator.py   # parse_clans_file(), generate_clans_file()
│   ├── similarity/
│   │   ├── tool_type.py              # ToolType enum (FOLDSEEK, USALIGN)
│   │   ├── struct_sim_tool.py        # Abstract base class; _execute_run()
│   │   ├── foldseek.py               # Foldseek: createdb → search → convertalis
│   │   ├── usalign.py                # USalign: direct subprocess call
│   │   └── struct_sim_computer.py    # StructSimComputer: factory + orchestration
│   ├── utils/
│   │   ├── fasta_utils.py            # FASTA parsing, UID extraction, region handling
│   │   ├── structure_utils.py        # AlphaFold download, region extraction (CIF)
│   │   ├── api_utils.py              # UniProt / UniParc REST calls
│   │   ├── file_utils.py             # reset_dir_content(), copy_dir_content(), download_file()
│   │   ├── dependency_checks.py      # verify_tool_dependencies()
│   │   └── log.py                    # setup_logging(verbose, quiet)
│   ├── benchmark/
│   │   ├── benchmark.py              # Benchmark class; run_all_tools(), export_csv()
│   │   └── benchmark_result.py       # BenchmarkResult dataclass
│   └── legacy/
│       ├── utils_old_clans.py        # headless CLANS Java execution (sequence-based)
│       ├── CLANS_website_version.jar
│       └── CLANS.jar
├── evaluation/                       # standalone analysis notebooks & modules
│   ├── evaluation_src/               # ScoresEvaluator, ClusterAnalyzer, ClansVisualizer, …
│   └── notebooks/
├── tests/
│   ├── conftest.py                   # shared fixtures (fixtures_dir, small_fasta_path, …)
│   ├── fixtures/                     # small.fasta, small.clans, small.tsv, small.cif
│   ├── unit/                         # pure unit tests; all I/O mocked
│   │   ├── core/
│   │   ├── similarity/
│   │   └── utils/
│   ├── integration/                  # pipeline tests; fixture files used; no real network/subprocess
│   ├── regression/                   # output shape stability tests
│   └── e2e/                          # require Foldseek/USalign in PATH + network access
└── .github/
    └── copilot-instructions.md       # GitHub Copilot-specific instructions (managed by GitHub; agents should not modify this file)
```

---

## 3. Environment Setup

### Prerequisites

- Python ≥ 3.10
- At least one external similarity tool on your `PATH`: `foldseek` and/or `USalign`

### Install (development mode)

```bash
git clone https://github.com/ElhabashyLab/3DClans.git
cd 3DClans

python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate

pip install -e ".[dev]"            # installs package + pytest, pytest-cov
```

> `pip install -e ".[dev]"` is required so that `import clans3d` resolves from `src/` and the `3dclans` CLI command is available.

### Verify the install

```bash
3dclans -h          # should print the CLI help

python -c "import clans3d; print('ok')"

which foldseek      # Linux / macOS
which USalign
```

---

## 4. Running Tests

All commands assume the venv is active and you are in the project root.

### Unit tests (fast, no external tools or network)

```bash
pytest tests/unit/
```

### Unit tests with coverage report

```bash
pytest tests/unit/ --cov=src/clans3d --cov-report=term-missing
```

### Integration tests (fixture files only, no real external calls)

```bash
pytest tests/integration/
```

### Regression tests (output shape stability)

```bash
pytest tests/regression/
```

### End-to-end tests (requires Foldseek/USalign + network)

```bash
pytest tests/e2e/ -m e2e
```

### Run a specific test file or class

```bash
pytest tests/unit/utils/test_fasta_utils.py
pytest tests/unit/utils/test_fasta_utils.py::TestExtractUidFromRecordID
```

### Test organization rules

| Test tier | Location | External I/O | Speed |
|-----------|----------|--------------|-------|
| Unit | `tests/unit/` | fully mocked | < 1 s per test |
| Integration | `tests/integration/` | fixture files; subprocess + network patched | seconds |
| Regression | `tests/regression/` | fixture files; no network | < 1 s per test |
| E2E | `tests/e2e/` | real network + real binaries | minutes |

Unit tests mirror the package layout: `tests/unit/core/`, `tests/unit/similarity/`, `tests/unit/utils/`. New unit test files must follow the same pattern and use `unittest.mock` (stdlib) for all I/O. The `benchmark/` and `legacy/` modules are excluded from coverage (see `pyproject.toml`) and do not have corresponding `tests/unit/` subdirectories — do not create them.

See [`docs/TESTING_GUIDE.md`](docs/TESTING_GUIDE.md) for a full table of what each existing test covers.

---

## 5. CLI Usage

```
3dclans -l <INPUT_FILE> -i <INPUT_TYPE> -t <TOOL> [OPTIONS]
```

| Argument | Values | Notes |
|----------|--------|-------|
| `-l / --load` | any path | path to the input file |
| `-i / --input_type` | `fasta` `a2m` `a3m` `tsv` | must match file extension |
| `-t / --tool` | `foldseek` `USalign` | case-sensitive |
| `-o / --out` | path | output file path or directory (default `output/clans_files/`) |
| `-s / --score` | `evalue` `TM` | Foldseek only; default `evalue` |
| `-c / --conf` | path | config file; CLI args override config |
| `-w / --workers` | int | parallel download threads (default 10) |
| `-v / --verbose` | flag | DEBUG logging |
| `-q / --quiet` | flag | ERROR-only logging; mutually exclusive with `-v` |

### Examples

```bash
# Foldseek with default E-value scoring
3dclans -l examples/small_fasta_files/5.fasta -i fasta -t foldseek

# Foldseek with TM-score
3dclans -l examples/small_fasta_files/5.fasta -i fasta -t foldseek -s TM

# USalign
3dclans -l examples/small_fasta_files/5.fasta -i fasta -t USalign

# TSV input
3dclans -l examples/small_tsv_files/5.tsv -i tsv -t foldseek

# Custom output path (file or directory)
3dclans -l examples/small_fasta_files/5.fasta -i fasta -t foldseek -o results/my_output.clans
3dclans -l examples/small_fasta_files/5.fasta -i fasta -t foldseek -o results/

# Config file
3dclans -c examples/config_files/example.conf
```

### Pipeline as a library

```python
from clans3d.core.pipeline import ClansPipeline, PipelineConfig
from clans3d.core.input_file_type import InputFileType
from clans3d.similarity.tool_type import ToolType

config = PipelineConfig(
    input_file="examples/small_fasta_files/5.fasta",
    input_type=InputFileType.FASTA,
    tool=ToolType.FOLDSEEK,
    output_dir="results",
    output_filename="my_output.clans",
)
pipeline = ClansPipeline(config)
clans_path, fasta_path = pipeline.run()
```

---

## 6. Architecture & Key Patterns

### Pipeline steps (`ClansPipeline`)

Each step is a separate public method, so benchmarking tools can time them independently:

| Step | Method | Description |
|------|--------|-------------|
| 1 | `fetch_structures()` | Download AlphaFold CIF files; returns `{uid: region | None}` |
| 2 | `generate_cleaned_fasta()` | Build cleaned FASTA for successfully downloaded structures |
| 3 | `compute_scores()` | Run similarity tool; returns DataFrame `[Sequence_ID_1, Sequence_ID_2, score]` |
| 4 | `generate_clans_file()` | Write `.clans` file; returns output path |

`pipeline.run()` chains all four steps.

### Enum-based type safety

Never compare tool or input type strings directly. Always use the enums:

```python
from clans3d.core.input_file_type import InputFileType
from clans3d.similarity.tool_type import ToolType

# Correct
if input_type == InputFileType.FASTA:
    ...

# Wrong — do not do this
if input_type == "fasta":
    ...
```

Convert from a string with `InputFileType("fasta")` or `ToolType("foldseek")`. An invalid string raises `ValueError`.

### Adding a similarity tool (StructSimTool subclass)

```python
# src/clans3d/similarity/mytool.py
import pandas as pd
from clans3d.similarity.struct_sim_tool import StructSimTool

class MyTool(StructSimTool):
    def __init__(self, working_dir: str):
        super().__init__("mytool", "My description", working_dir)

    def start_run(self, structures_dir: str, expected_number_of_scores: int) -> pd.DataFrame:
        self.command = ["mytool", "--input", structures_dir, "--output", self.working_dir]
        return self._execute_run(expected_number_of_scores)

    def _log_progress(self, stdout_reader: dict, stderr_reader: dict) -> None:
        # Log a progress message using stdout_reader["latest_line"] or stderr_reader["latest_line"]
        pass

    def _parse_output(self) -> pd.DataFrame:
        # Parse self.output (captured stdout string) into a DataFrame
        # Must return columns: [Sequence_ID_1, Sequence_ID_2, score]
        ...
```

Register the new tool in `StructSimComputer._create_tool()` and add a value to `ToolType`.

### CLANS file sections

```
sequences=N
<param>
key=value
</param>
<seq>
>header
SEQUENCE
</seq>
<pos>
idx x y z
</pos>
<hsp>
idx1 idx2:score
</hsp>
```

- `<pos>` and `<hsp>` use **integer indices** that correspond to the order of sequences in `<seq>`.
- Internally (`ClansFile` and `ClansFileGenerator`) UIDs are used throughout; conversion to/from indices happens automatically. See the [Common Gotchas](#11-common-gotchas) section for details.

### Structure files

- Format: **CIF (mmCIF)** only — no PDB format used.
- Stored in `work/structures/<uid>.cif`.
- Region extraction: `extract_region_of_protein()` in `structure_utils.py` uses BioPython `MMCIFParser` + `MMCIFIO`.
- Both Foldseek and USalign accept CIF natively; no conversion needed.

### Directory layout at runtime

```
work/
├── structures/                  # downloaded .cif files (reset per run)
├── cleaned_input_storage/       # <input_name>_cleaned.fasta
├── foldseek/                    # Foldseek working files (reset per run)
└── usalign/                     # USalign working files (reset per run)
output/
└── clans_files/
    └── <input_name>_cleaned.clans
```

`reset_dir_content(dir)` deletes and recreates a directory. Tool working directories are always reset at the start of a run.

---

## 7. Data Contracts

### Pairwise scores DataFrame

All similarity tools must return a DataFrame with exactly these columns:

| Column | Type | Description |
|--------|------|-------------|
| `Sequence_ID_1` | `str` | UniProt accession (UID) of the first protein |
| `Sequence_ID_2` | `str` | UniProt accession (UID) of the second protein |
| `score` | `float` | Similarity score (e-value or `1 - max(TM1, TM2)`) |

- Self-hits (`Sequence_ID_1 == Sequence_ID_2`) must be removed.
- Symmetric duplicates (A→B and B→A) must be deduplicated — keep one row per pair.
- Expected count: `n * (n - 1) / 2` where `n` is the number of structures.

### `fetch_structures` return value

```python
dict[str, tuple[int, int] | None]
# uid → (region_start, region_end) if a region was specified, else None
```

### FASTA record ID conventions

- Full header format: `sp|ACCESSION|ENTRY_NAME/start-end`
- `extract_uid_from_recordID(record_id)` returns the bare UniProt accession (e.g. `P49811`).
- `extract_region_from_record(record)` returns `(start, end)` or `None`.
- Region notation is 1-based, inclusive.

---

## 8. Adding New Features

### New similarity tool

1. Create `src/clans3d/similarity/<toolname>.py` — subclass `StructSimTool`.
2. Implement `start_run()`, `_log_progress()`, and `_parse_output()`.
3. Add a value to `ToolType` enum in `src/clans3d/similarity/tool_type.py`.
4. Register the tool in `StructSimComputer._create_tool()`.
5. Add unit tests in `tests/unit/similarity/test_<toolname>.py` — mock the subprocess.
6. Update `src/clans3d/utils/dependency_checks.py` if the tool needs a binary check.

### New input format

1. Add a value to `InputFileType` in `src/clans3d/core/input_file_type.py`.
2. Handle the new type in `fetch_structures()` in `structure_utils.py`.
3. Handle cleaned FASTA generation in `ClansPipeline.generate_cleaned_fasta()`.
4. Add the new extension to `ClansPipeline._validate()`.
5. Update `cli.py` `choices` list for `--input_type`.
6. Add unit tests.

### New config option

1. Add the argument to `_build_main_parser()` in `core/cli.py`.
2. Add the corresponding attribute to `PipelineConfig` in `core/pipeline.py`.
3. Pass the value through `main()` in `main.py`.
4. Update documentation (see [Documentation updates](#documentation-updates) below).

### Documentation updates

**Every user-facing change must be reflected in the documentation.** After implementing any feature, bug fix, or behavioral change, check and update the following files as needed:

- **`README.md`** — CLI argument tables, usage examples, output descriptions, configuration file examples.
- **`AGENTS.md`** — CLI usage table (section 5), examples, pipeline-as-a-library snippets, "Adding New Features" checklists.
- **`docs/TESTING_GUIDE.md`** — if new test files or test classes were added.
- **`docs/BENCHMARK_USAGE.md`** — if benchmark behavior or output changed.

When in doubt, update the docs. Outdated documentation is worse than missing documentation.

---

## 9. Code Conventions

### Naming

- All Python source files: `snake_case.py`
- Classes: `PascalCase`
- Functions/variables: `snake_case`
- Enum files live beside the module they belong to: `tool_type.py` in `similarity/`, `input_file_type.py` in `core/`

### Imports

Always use full `clans3d` package paths — never relative imports:

```python
# Correct
from clans3d.core.input_file_type import InputFileType

# Wrong
from .input_file_type import InputFileType
```

### Logging

Use `logging.getLogger(__name__)` at module level. Never use `print()` for status messages in library code. Use the appropriate level:

- `logger.debug()` — detailed internals (shown only with `-v`)
- `logger.info()` — normal pipeline progress
- `logger.warning()` — recoverable issues (e.g., download failure for one entry)
- `logger.error()` — non-recoverable errors (caught in `main()`)

### Error handling

- `FileNotFoundError` — missing input file
- `ValueError` — invalid argument, missing required columns, malformed data
- `RuntimeError` — subprocess failure, no structures downloaded
- Warn on score count mismatches (do not fail); the `StructSimComputer` logs a warning.

### Mocking in tests

Use `unittest.mock` (stdlib). Patch at the call site, not at the definition site:

```python
# Patch where it is used, not where it is defined
with patch("clans3d.core.pipeline.fetch_structures") as mock_fetch:
    mock_fetch.return_value = {"P11111": None}
    ...
```

External calls that must always be mocked in unit/integration tests:
- `requests.get` / any network I/O
- `subprocess.run` / `subprocess.Popen`
- `download_file()` from `file_utils`
- `fetch_structures()` from `structure_utils`

### Comments

Only add comments when they match the style of the surrounding code or when explaining non-obvious logic. Do not add docstrings to trivial one-liners.

---

## 10. External Dependencies

### Python packages (installed via `pip install -e ".[dev]"`)

| Package | Purpose |
|---------|---------|
| `biopython` | FASTA parsing, CIF parsing, SeqRecord |
| `pandas` | Pairwise score DataFrames |
| `requests` | AlphaFold API, UniProt REST calls |
| `scipy`, `scikit-learn`, `hdbscan` | Evaluation clustering |
| `networkx`, `python-igraph`, `leidenalg` | Graph-based clustering |
| `matplotlib`, `seaborn` | Visualization |
| `pytest`, `pytest-cov` | Testing (dev only) |

### External binaries (must be on PATH)

| Binary | Install | Notes |
|--------|---------|-------|
| `foldseek` | https://mmseqs.com/foldseek | Linux / macOS only |
| `USalign` | https://zhanggroup.org/US-align/ | Linux / macOS / Windows |

`verify_tool_dependencies(tool_type)` in `dependency_checks.py` checks for the binary at startup and raises `RuntimeError` if it is missing.

---

## 11. Common Gotchas

### UID extraction is non-trivial

`extract_uid_from_recordID()` uses a regex against the official UniProt accession pattern. It raises `ValueError` for any string that does not contain a valid accession. Test with both `sp|P49811|...` and plain `P49811` headers.

### Foldseek score types behave differently

- `evalue`: lower is more similar (closer to 0 → stronger hit); used directly as the score.
- `TM`: two TM-scores are returned (`qtmscore`, `ttmscore`); the final score is `1 - max(TM1, TM2)` — **lower means more similar**.

Both Foldseek score types and USalign output `score` column values where **lower = more similar**. Keep this in mind when writing threshold logic.

### Input type must match file extension

`ClansPipeline._validate()` checks that the `-i` flag matches the file extension. Mismatches raise `ValueError` immediately. The extension check is case-insensitive.

### Working directories are reset on every run

`reset_dir_content(dir)` deletes and recreates the directory. Never store files in `work/` that you need to persist across pipeline runs.

### Config file format uses single dashes

Config files use `-key value` (one dash), but the equivalent CLI argument uses `--key value` (two dashes). `config_to_argv()` converts automatically.

### E2E tests are opt-in

E2E tests are marked `@pytest.mark.e2e` and are skipped automatically when Foldseek/USalign are not installed. To run them you must explicitly pass `-m e2e`:

```bash
pytest tests/e2e/ -m e2e
```

Never add E2E-style network or subprocess calls to unit or integration tests.

### Benchmark and legacy modules are excluded from coverage

The coverage configuration in `pyproject.toml` explicitly omits `benchmark/`, `legacy/`, and `dataset_generator/`. Do not add coverage expectations for those modules.

### ClansFile uses UIDs internally; indices only in the serialized file

`ClansFile.__str__()` converts UIDs → integer indices automatically via `_uid_to_index`. Do not manually insert integer indices into `coordinates` or `scores` DataFrames that are passed to `ClansFile`.
