# Testing Guide for Clans-3D

## Overview

This document outlines the testing strategy for Clans-3D. Because the application integrates external binaries (Foldseek, USalign), network APIs (UniProt), and file I/O-heavy workflows, tests are split into several layers with different requirements and run frequencies.

---

## Test Types

### 1. Unit Tests

**What:** Test individual functions and classes in isolation, with all external dependencies mocked or replaced with fixtures.

**Tooling:** `pytest` + `unittest.mock`

**What to cover:**

| Module | Test targets |
|---|---|
| `utils/fasta_utils.py` | `extract_uid_from_recordID()`, `extract_uids_from_fasta()`, `clean_fasta_file()`, region parsing |
| `core/clans_file_generator.py` | `parse_clans_file()`, `generate_clans_file_*()` — round-trip test on a known `.clans` fixture |
| `core/clans_file.py` | ClansFile dataclass construction and field access |
| `core/config_file.py` | Config parse/write round-trip |
| `evaluation/data_normalizer.py` | min-max, z-score, -log10 on known DataFrames with expected outputs |
| `evaluation/cluster_analyzer.py` | DBSCAN/HDBSCAN/Leiden on small synthetic coordinate DataFrames |
| `similarity/struct_sim_computer.py` | `_create_tool()` factory — correct subclass returned per `ToolType` |
| `similarity/struct_sim_tool.py` | Score-count mismatch warning logic |
| `core/input_file_type.py` | Enum value mapping |

**Example:**

```python
# tests/unit/test_fasta_utils.py
import pytest
from clans3d.utils.fasta_utils import extract_uid_from_recordID

@pytest.mark.parametrize("record_id, expected", [
    ("sp|P49811|MYOD1_PIG", "P49811"),
    ("sp|P49811|MYOD1_PIG/2-300", "P49811"),
    ("P49811", "P49811"),
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
├── conftest.py                  # Shared fixtures (tmp_path wrappers, small DataFrames, etc.)
├── fixtures/
│   ├── small.fasta              # 5-entry FASTA for fast tests
│   ├── small.clans              # Pre-generated CLANS file matching small.fasta
│   └── scores_5.tsv             # Pre-computed pairwise scores for 5 sequences
├── unit/
│   ├── test_fasta_utils.py
│   ├── test_clans_file_generator.py
│   ├── test_data_normalizer.py
│   ├── test_cluster_analyzer.py
│   └── test_struct_sim_computer.py
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

1. **Unit tests** for `fasta_utils`, `clans_file_generator`, and `data_normalizer` — pure functions with no external dependencies; highest ROI
2. **Integration tests** for the CLANS round-trip and pipeline with mocked tool execution
3. **Regression tests** to lock in output stability before refactoring
4. **Contract tests** for API calls to decouple from network availability
5. **E2E tests** for final validation with real tools
6. **Property-based tests** for normalizer invariants and UID extraction robustness

---

## Tooling Summary

| Library | Purpose | Install |
|---|---|---|
| `pytest` | Test runner | `pip install pytest` |
| `pytest-cov` | Coverage reporting | `pip install pytest-cov` |
| `unittest.mock` | Mocking (stdlib) | built-in |
| `responses` | HTTP mocking | `pip install responses` |
| `syrupy` | Snapshot/regression testing | `pip install syrupy` |
| `hypothesis` | Property-based testing | `pip install hypothesis` |

Run all non-E2E tests:

```bash
pytest tests/ -m "not e2e" --cov=src/clans3d --cov-report=term-missing
```

Run only E2E tests (requires Foldseek/USalign in PATH):

```bash
pytest tests/e2e/ -m e2e
```
