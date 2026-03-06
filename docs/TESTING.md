# Testing Guide for Clans-3D

## Running Tests

All test commands below assume you are in the project root and the virtual environment is active.

### Unit tests (fast, no external tools or network required)

```bash
pytest tests/unit/
```

With coverage report:

```bash
pytest tests/unit/ --cov=src/clans3d --cov-report=term-missing
```

Target a specific module:

```bash
pytest tests/unit/utils/test_fasta_utils.py
```

Run a single test by name:

```bash
pytest tests/unit/utils/test_fasta_utils.py::test_extract_uid_from_record_id
```

### All non-E2E tests (unit + integration + regression)

```bash
pytest tests/ -m "not e2e"
```

### E2E tests (requires Foldseek / USalign in PATH)

```bash
pytest tests/ -m e2e
```

> **Note:** E2E tests are not yet implemented. The `-m e2e` flag and the `e2e/` directory are reserved for them.

---

## Overview

This document outlines the testing strategy for Clans-3D. Because the application integrates external binaries (Foldseek, USalign), network APIs (UniProt), and file I/O-heavy workflows, tests are split into several layers with different requirements and run frequencies.

---

## Test Types

### 1. Unit Tests

**What:** Test individual functions and classes in isolation, with all external dependencies mocked or replaced with fixtures.

**Tooling:** `pytest` + `unittest.mock`

**What to cover:**

#### `utils/fasta_utils.py`

| Function                       | What to test                                                                                                                                             |
| ------------------------------ | -------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `extract_uid_from_recordID()`  | Standard `sp\|UID\|...` format, `tr\|UID\|...`, plain UID, UID with version suffix (`P12345.1`), UID with region `/1-68`, underscore-style `UID_SPECIES` |
| `extract_region_from_record()` | Header with region (`/1-68`), header without region, region with large numbers                                                                           |
| `add_region_to_record()`       | With valid region trims sequence and updates ID; `None` region leaves record unchanged; invalid region (out of bounds, start > end) raises `ValueError`  |
| `create_mock_up_record()`      | Correct ID format with and without region; sequence is `"not_found"`                                                                                     |
| `extract_uids_from_fasta()`    | Reads a small fixture FASTA with known UIDs; returns them in order                                                                                       |
| `extract_records_from_fasta()` | Returns correct number of `SeqRecord` objects with UIDs as IDs                                                                                           |
| `clean_fasta_file()`           | Retains only records whose UIDs are in the allow-list; writes to file                                                                                    |
| `copy_records_from_fasta()`    | Only records matching given UIDs are written                                                                                                             |

#### `utils/api_utils.py`

| Function                                     | What to test                                                                                               |
| -------------------------------------------- | ---------------------------------------------------------------------------------------------------------- |
| `_chunked()`                                 | Even split, uneven split, empty list, chunk size larger than list                                          |
| `uniprot_accessions_to_uniparc_accessions()` | Mock all HTTP calls; assert correct UID→UPI mapping returned; assert `None` for UIDs not found in response |

#### `utils/dependency_checks.py`

| Function                     | What to test                                                                                                   |
| ---------------------------- | -------------------------------------------------------------------------------------------------------------- |
| `check_external_tool()`      | Returns a path string when tool exists (`shutil.which` mocked to return a path); returns `None` when not found |
| `verify_tool_dependencies()` | Raises `RuntimeError` when tool not in PATH; does not raise when tool is present                               |

#### `utils/file_utils.py`

| Function              | What to test                                                                                                     |
| --------------------- | ---------------------------------------------------------------------------------------------------------------- |
| `reset_dir_content()` | Directory created if missing; existing files removed before recreation                                           |
| `copy_dir_content()`  | Files and nested subdirectories are copied; target is created if missing                                         |
| `download_file()`     | Returns `True` and writes file on HTTP 200 (mock `requests.get`); returns `False` and logs warning on HTTP error |

#### `utils/structure_utils.py`

| Function                         | What to test                                                                                                                                                           |
| -------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `_log_interval()`                | Returns correct interval for small, medium, and large totals                                                                                                           |
| `extract_region_of_protein()`    | Given a real or minimal mock CIF file + region, output contains only residues within range                                                                             |
| `download_alphafold_structure()` | With mocked API: returns `True` when URL resolves and file downloads; returns `False` on API error; calls `extract_region_of_protein()` only when region is not `None` |
| `fetch_structures()`             | Dispatches to `process_fasta_file` for FASTA/A2M input types; dispatches to `process_tsv_file` for TSV; raises `ValueError` for unsupported type                       |

#### `utils/log.py`

| Function          | What to test                                                                                                                  |
| ----------------- | ----------------------------------------------------------------------------------------------------------------------------- |
| `setup_logging()` | `verbose=True` → `DEBUG` level; `quiet=True` → `ERROR` level; default → `INFO`; calling twice does not add duplicate handlers |

---

#### `core/clans_file.py`

| Aspect                         | What to test                                                                                                                          |
| ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------- |
| `__init__()`                   | Raises `ValueError` when `fasta_records` and `path_to_fasta` are both `None`; accepts either                                          |
| `_validate_data_consistency()` | Raises when coordinate count ≠ sequence count; raises when coordinate UIDs not in records; raises when score UIDs not in records      |
| `__str__()`                    | Output contains `sequences=N`, `<param>`, `<seq>`, `<pos>`, `<hsp>` sections; indices in `<pos>`/`<hsp>` match the FASTA record order |
| `_uid_to_index` mapping        | Built correctly from fasta_records in order                                                                                           |

#### `core/clans_file_generator.py`

| Method                           | What to test                                                                        |
| -------------------------------- | ----------------------------------------------------------------------------------- |
| `_extract_block()`               | Returns lines between tags; returns `None` when tag missing                         |
| `_parse_number_of_sequences()`   | Correct parsing; raises `ValueError` when `sequences=` line is absent               |
| `_parse_param_block()`           | Returns dict of params; handles absent param block (returns `None`)                 |
| `_parse_fasta_block()`           | Returns correct number of `SeqRecord` objects; raises on missing `<seq>` block      |
| `_parse_pos_block()`             | Returns correct (idx, x, y, z) tuples; raises on missing/malformed lines            |
| `_parse_scores_block()`          | Returns DataFrame with correct columns and dtypes; handles empty `<hsp>` block      |
| `_normalize_scores_format()`     | Removes duplicate A↔B / B↔A pairs; result has fewer rows than input with duplicates |
| `_generate_random_coordinates()` | Returns one tuple per UID; x/y/z all in [0, 1]                                      |

#### `core/config_file.py`

| Method             | What to test                                                                                                  |
| ------------------ | ------------------------------------------------------------------------------------------------------------- |
| `write_config()`   | Written file contains `-key value` lines for each entry; documentation line written as comment                |
| `read_config()`    | Round-trips a written config; blank lines and `#` comments ignored; raises `ValueError` on missing `-` prefix |
| `config_to_argv()` | Returns `["--key", "value", ...]` list matching config contents                                               |

#### `core/pipeline.py`

| Aspect                      | What to test                                                                                                                                                      |
| --------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `ClansPipeline._validate()` | Raises `FileNotFoundError` for missing input file; raises `ValueError` on extension mismatch; raises `ValueError` when `foldseek_score` set for non-Foldseek tool |
| `PipelineConfig` defaults   | `structures_dir`, `output_dir`, `cleaned_input_storage` use expected default paths                                                                                |

#### `core/cli.py`

| Aspect         | What to test                                                                                                                                              |
| -------------- | --------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `parse_args()` | All required args present → namespace with correct values; config file values merged and overridden by CLI args; invalid tool or input type → parse error |

---

#### `similarity/tool_type.py` and `core/input_file_type.py`

| Aspect          | What to test                                                                                 |
| --------------- | -------------------------------------------------------------------------------------------- |
| Enum round-trip | `ToolType("foldseek") == ToolType.FOLDSEEK`; `InputFileType("fasta") == InputFileType.FASTA` |
| Invalid value   | `ToolType("invalid")` raises `ValueError`                                                    |

#### `similarity/struct_sim_computer.py`

| Method           | What to test                                                                                                       |
| ---------------- | ------------------------------------------------------------------------------------------------------------------ |
| `_create_tool()` | Returns `Foldseek` for `ToolType.FOLDSEEK`; `USalign` for `ToolType.USALIGN`; raises `ValueError` for unknown type |
| `run()`          | Logs warning when returned score count ≠ expected (patch `tool.start_run` to return short DataFrame)               |

#### `similarity/foldseek.py`

| Method                                  | What to test                                                                                                         |
| --------------------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| `_get_output_columns_from_self_score()` | Returns `"query,target,evalue"` for `"evalue"`; returns TM columns for `"TM"`; raises `ValueError` for unknown score |
| `_clean_scores()` with `score="evalue"` | Self-hits removed; symmetric duplicates removed; `score` column equals original `evalue`                             |
| `_clean_scores()` with `score="TM"`     | `score = 1 - max(TM1, TM2)`; output has `score` column, no `TM1`/`TM2`                                               |
| `_detect_phase()`                       | Returns correct label when known keyword present in chunks; returns `""` when no match                               |

#### `similarity/usalign.py`

| Method                | What to test                                                                                                                                        |
| --------------------- | --------------------------------------------------------------------------------------------------------------------------------------------------- |
| `_count_data_lines()` | Ignores lines starting with `#` or `Warning`; counts remaining non-empty lines                                                                      |
| `_parse_output()`     | Mock `self.output` with known TSV string; assert resulting DataFrame has correct columns, `score = 1 - max(TM1, TM2)`, extensions stripped from IDs |

#### `evaluation/data_normalizer.py`

| Method                | What to test                                                                         |
| --------------------- | ------------------------------------------------------------------------------------ |
| min-max normalization | Output min is 0.0, output max is 1.0 for non-constant input; constant input handling |
| z-score normalization | Output mean ≈ 0, std ≈ 1 for known input                                             |
| -log10 normalization  | Correct transform for known e-value inputs; handles 0 or negative values gracefully  |

#### `evaluation/cluster_analyzer.py`

| Method                          | What to test                                                                                                                |
| ------------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
| `find_clusters_density_based()` | On small synthetic 3D coordinates with obvious clusters: DBSCAN assigns correct number of clusters; noise points labeled -1 |
| `find_clusters_graph_based()`   | Leiden on a simple synthetic score DataFrame assigns expected community structure                                           |

**Example:**

```python
# tests/unit/test_fasta_utils.py
import pytest
from clans3d.utils.fasta_utils import extract_uid_from_recordID

@pytest.mark.parametrize("record_id, expected", [
    ("sp|P49811|MYOD1_PIG", "P49811"),
    ("sp|P49811|MYOD1_PIG/2-300", "P49811"),
    ("tr|A0A2M8UWB6.1|A0A2M8UWB6_PSESP/524-595", "A0A2M8UWB6"),
    ("P49811", "P49811"),
    ("P49811.1", "P49811"),
    ("MYOD1_PIG", "MYOD1"),
])
def test_extract_uid_from_record_id(record_id, expected):
    assert extract_uid_from_recordID(record_id) == expected
```

---

### 2. Integration Tests

**What:** Test that multiple components work together correctly end-to-end, using real files from `tests/fixtures/` but without executing external binaries or making network calls.

**Tooling:** `pytest` + `tmp_path` fixture + pre-computed score fixtures

**What to cover:**

- `ClansFileGenerator.parse_clans_file()` → `ClansFileGenerator.generate_clans_file_*()` round-trip using the real `.clans` files in `examples/clans_files/`
- `StructSimComputer.run()` with a mocked `StructSimTool` that returns a known DataFrame, verifying the score-count check and return value
- `ScoresEvaluator.extract_data_from_clans_files()` merging structural and sequence score DataFrames correctly
- Full pipeline (`ClansPipeline.run()`) using patched `fetch_structures()` and patched tool execution returning pre-computed scores

**Example:**

```python
# tests/integration/test_clans_round_trip.py
from clans3d.core.clans_file_generator import ClansFileGenerator

def test_parse_generate_round_trip(tmp_path):
    generator = ClansFileGenerator()
    clans = generator.parse_clans_file("examples/clans_files/500_cleaned.clans")
    out = tmp_path / "out.clans"
    generator.generate_clans_file_from_clans(clans, str(out))
    clans2 = generator.parse_clans_file(str(out))
    assert len(clans.fasta_records) == len(clans2.fasta_records)
    assert clans.scores_df.shape == clans2.scores_df.shape
```

---

### 3. End-to-End (E2E) / System Tests

**What:** Execute the real CLI or `ClansPipeline.run()` against a small input file with actual external tools. These are slow and require Foldseek and/or USalign to be installed.

**Tooling:** `pytest` with a `--run-e2e` flag or a separate `pytest` mark (`@pytest.mark.e2e`)

**What to cover:**

- `clans3d -l examples/small_fasta_files/5.fasta -i fasta -t foldseek` produces a valid `.clans` output file
- Output file parses without error and has the expected number of sequences and score entries
- Same scenario with USalign

**Tip:** Gate these tests behind a mark so they don't run in CI unless explicitly requested:

```python
@pytest.mark.e2e
def test_full_pipeline_foldseek(tmp_path):
    ...
```

```ini
# pytest.ini or pyproject.toml
[pytest]
markers =
    e2e: end-to-end tests requiring external tools
```

---

### 4. Regression Tests

**What:** Lock in known-good outputs so that refactors don't silently change behavior.

**Tooling:** `pytest` + snapshot testing via `syrupy` or manual baseline file comparison

**What to cover:**

- For a fixed small input (e.g. `examples/small_fasta_files/5.fasta`) with pre-saved mock scores, the generated `.clans` file matches a stored baseline byte-for-byte (or structurally via parsed comparison)
- `DataNormalizer` output values for a fixed DataFrame do not drift
- Cluster labels for a fixed set of 2D/3D coordinates remain stable across algorithm versions

---

### 5. Contract / API Tests

**What:** Verify that the external HTTP calls to UniProt/UniParc/AlphaFold return the expected shape of data (not necessarily the exact data). Use `responses` or `pytest-httpserver` to mock the HTTP layer.

**Tooling:** `pytest` + `responses` library

**What to cover:**

- `api_utils.uniprot_accessions_to_uniparc_accessions()` — mock the REST endpoint, assert correct column extraction
- `structure_utils.fetch_structures()` — mock the AlphaFold CIF download, assert files are written to the expected paths with the expected names

---

### 6. Property-Based Tests

**What:** Generate random but valid inputs and assert that invariants always hold, catching edge cases not covered by handwritten examples.

**Tooling:** `hypothesis`

**Good candidates:**

- `extract_uid_from_recordID()` — any string of the form `sp|XXXXX|...` should return the second pipe-delimited segment
- `DataNormalizer` — min-max output is always in `[0, 1]` for any non-constant input
- Score DataFrame merging — merged output always has `≤ min(len(df1), len(df2))` rows when doing an inner join

---

## Recommended Directory Layout

```
tests/
├── __init__.py
├── conftest.py                       # Shared fixtures (tmp_path wrappers, small DataFrames, etc.)
├── fixtures/
│   ├── small.fasta                   # 5-entry FASTA for fast tests
│   ├── small.clans                   # Pre-generated CLANS file matching small.fasta
│   ├── scores_5.tsv                  # Pre-computed pairwise scores for 5 sequences
│   └── small.cif                     # Minimal CIF structure for region extraction tests
├── unit/
│   ├── utils/
│   │   ├── test_fasta_utils.py
│   │   ├── test_api_utils.py
│   │   ├── test_dependency_checks.py
│   │   ├── test_file_utils.py
│   │   ├── test_structure_utils.py
│   │   └── test_log.py
│   ├── core/
│   │   ├── test_clans_file.py
│   │   ├── test_clans_file_generator.py
│   │   ├── test_config_file.py
│   │   ├── test_pipeline_validate.py
│   │   ├── test_cli.py
│   │   └── test_input_file_type.py
│   ├── similarity/
│   │   ├── test_tool_type.py
│   │   ├── test_struct_sim_computer.py
│   │   ├── test_foldseek.py
│   │   └── test_usalign.py
│   └── evaluation/
│       ├── test_data_normalizer.py
│       └── test_cluster_analyzer.py
├── integration/
│   ├── test_clans_round_trip.py
│   ├── test_pipeline_mocked.py
│   └── test_scores_evaluator.py
├── e2e/
│   ├── test_pipeline_foldseek.py
│   └── test_pipeline_usalign.py
└── regression/
    └── test_output_stability.py
```

---

## Suggested Priority Order

1. **Unit tests — pure utility functions** (`fasta_utils`, `file_utils`, `config_file`, `data_normalizer`, enum types) — no I/O or subprocess; fastest to write and run
2. **Unit tests — parsing logic** (`clans_file_generator` private parsers, `ClansFile.__str__` / `_validate_data_consistency`, `USalign._parse_output`, `Foldseek._clean_scores`) — self-contained string/DataFrame logic
3. **Unit tests — mocked I/O** (`dependency_checks`, `structure_utils`, `api_utils`, `ClansPipeline._validate`) — require `unittest.mock` for `shutil.which`, `requests`, `os.path`
4. **Integration tests** for the CLANS round-trip (`parse_clans_file` → `generate_clans_file` with real fixture files) and pipeline with mocked tool execution
5. **Regression tests** to lock in output stability before refactoring
6. **Contract tests** for API calls to decouple from network availability
7. **E2E tests** for final validation with real tools installed
8. **Property-based tests** for normalizer invariants and UID extraction robustness

---

## Tooling Summary

| Library         | Purpose                     | Install                  |
| --------------- | --------------------------- | ------------------------ |
| `pytest`        | Test runner                 | `pip install pytest`     |
| `pytest-cov`    | Coverage reporting          | `pip install pytest-cov` |
| `unittest.mock` | Mocking (stdlib)            | built-in                 |
| `responses`     | HTTP mocking                | `pip install responses`  |
| `syrupy`        | Snapshot/regression testing | `pip install syrupy`     |
| `hypothesis`    | Property-based testing      | `pip install hypothesis` |

Run all non-E2E tests:

```bash
pytest tests/ -m "not e2e" --cov=src/clans3d --cov-report=term-missing
```

Run only E2E tests (requires Foldseek/USalign in PATH):

```bash
pytest tests/e2e/ -m e2e
```
