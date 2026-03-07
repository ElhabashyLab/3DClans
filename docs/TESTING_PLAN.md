# Testing Plan for Clans-3D

This document describes the additional test layers implemented beyond unit tests, and what each covers. Evaluation and benchmark code, as well as CLANS file parsing (used only by the evaluation layer), are out of scope.

---

## 1. Integration Tests

**Goal:** Verify that pipeline components work together end-to-end using real fixture files, without calling external binaries or making network calls.

**Location:** `tests/integration/`

**Tooling:** `pytest` + `tmp_path` fixture + `unittest.mock`

### Full pipeline with mocked structures and scores

Patch `fetch_structures()` and `StructSimTool.start_run()` to avoid any I/O or subprocess calls. Run `ClansPipeline.run()` end-to-end with `tests/fixtures/small.fasta`. Assert that:

- An output `.clans` file is written
- The file contains the expected number of sequences and score entries

### Run integration tests

```bash
pytest tests/integration/
```

---

## 2. Regression Tests

**Goal:** Lock in known-good structural properties of generated CLANS files so refactors cannot silently break output shape.

**Location:** `tests/regression/`

**Tooling:** `pytest` + `tmp_path` + mock scores from fixtures

### What is locked in

| Test                          | What is verified                                                                                                                                                        |
| ----------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `test_clans_output_structure` | For `tests/fixtures/small.fasta` with pre-built mock scores, the generated `.clans` file always has the correct sequence count, coordinate count, and score entry count |

### Run regression tests

```bash
pytest tests/regression/
```

---

## 3. End-to-End (E2E) Tests

**Goal:** Execute `ClansPipeline.run()` and the CLI against a small input file using the real Foldseek and USalign binaries. These tests are slow and are skipped automatically when the tools are not in `PATH`.

**Location:** `tests/e2e/`

**Tooling:** `pytest` + `@pytest.mark.e2e` + `shutil.which` skip guards

### Scenarios covered

| Test                         | Tool                | Input                                | What is asserted                                                            |
| ---------------------------- | ------------------- | ------------------------------------ | --------------------------------------------------------------------------- |
| `TestPipelineFoldseekEvalue` | Foldseek (evalue)   | `examples/small_fasta_files/5.fasta` | `.clans` file written; ≥ 3 sequences; ≥ 1 score entry; all sections present |
| `TestPipelineFoldseekTM`     | Foldseek (TM score) | `examples/small_fasta_files/5.fasta` | File written; all scores in `[0, 1]`                                        |
| `TestPipelineUSalign`        | USalign             | `examples/small_fasta_files/5.fasta` | File written; ≥ 3 sequences; ≥ 1 score entry                                |
| `TestCliEntrypoint`          | Foldseek            | `examples/small_fasta_files/5.fasta` | `clans3d` subprocess exits with code `0`                                    |

### Run E2E tests

```bash
pytest tests/e2e/ -m e2e
```

---

## Directory Layout

```
tests/
├── __init__.py
├── conftest.py                 # Shared fixtures + e2e skip guards
├── fixtures/
│   ├── small.fasta
│   ├── small.clans
│   ├── small.tsv
│   ├── small.cif
│   └── small_roi.cif
├── unit/                       # Implemented — see TESTING_GUIDE.md
│   ├── utils/
│   ├── core/
│   └── similarity/
├── integration/
│   └── test_pipeline_mocked.py
├── regression/
│   └── test_clans_output.py
└── e2e/
    └── test_pipeline_e2e.py
```

---

## Tooling Summary

| Library         | Purpose            | Install                  |
| --------------- | ------------------ | ------------------------ |
| `pytest`        | Test runner        | `pip install pytest`     |
| `pytest-cov`    | Coverage reporting | `pip install pytest-cov` |
| `unittest.mock` | Mocking (stdlib)   | built-in                 |

Run all non-E2E tests:

```bash
pytest tests/ -m "not e2e" --cov=src/clans3d --cov-report=term-missing
```

**Gate condition:** Skip automatically unless `foldseek` / `usalign` is found in `PATH`.

```python
# tests/conftest.py (addition)
import shutil, pytest

def pytest_configure(config):
    config.addinivalue_line("markers", "e2e: end-to-end tests requiring external tools")

@pytest.fixture
def require_foldseek():
    if not shutil.which("foldseek"):
        pytest.skip("foldseek not in PATH")

@pytest.fixture
def require_usalign():
    if not shutil.which("USalign"):
        pytest.skip("USalign not in PATH")
```

### Scenarios to cover

| Test                            | Tool                | Input                                | What to assert                                                           |
| ------------------------------- | ------------------- | ------------------------------------ | ------------------------------------------------------------------------ |
| `test_pipeline_foldseek_evalue` | Foldseek            | `examples/small_fasta_files/5.fasta` | Output `.clans` file exists, parses, has 5 sequences and ≥ 1 score entry |
| `test_pipeline_foldseek_TM`     | Foldseek (TM score) | `examples/small_fasta_files/5.fasta` | Same; scores are in `[0, 1]`                                             |
| `test_pipeline_usalign`         | USalign             | `examples/small_fasta_files/5.fasta` | Same as above                                                            |
| `test_cli_entrypoint`           | Foldseek            | `examples/small_fasta_files/5.fasta` | `subprocess.run(["clans3d", ...])` exits with code `0`                   |

Run only E2E tests:

```bash
pytest tests/e2e/ -m e2e
```

---

## 4. Regression Tests

**Goal:** Lock in known-good outputs so that refactors do not silently change behavior.

**Location:** `tests/regression/`

**Tooling:** `pytest` + file comparison or `syrupy` snapshot testing

### Scenarios to cover

| Test                           | What to lock in                                                                                                                                                                                           |
| ------------------------------ | --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| CLANS file output stability    | For a fixed input (e.g. `tests/fixtures/small.fasta`) with pre-saved mock scores, the generated `.clans` file matches a stored baseline structurally (parsed record count, score count, coordinate count) |
| `DataNormalizer` output values | For a fixed DataFrame input, normalized output values do not drift across algorithm changes                                                                                                               |
| Cluster label stability        | For a fixed set of 3D coordinates, DBSCAN/HDBSCAN/Leiden cluster labels remain stable across library version upgrades                                                                                     |

When using `syrupy`:

```bash
pip install syrupy
pytest tests/regression/ --snapshot-update   # update baselines
pytest tests/regression/                      # assert no drift
```

---

## 5. Contract / HTTP API Tests

**Goal:** Verify that the external HTTP calls to UniProt/UniParc/AlphaFold return data in the expected shape, decoupling tests from live network availability.

**Location:** `tests/unit/utils/` or a dedicated `tests/contract/`

**Tooling:** `pytest` + `responses` library

```bash
pip install responses
```

### Scenarios to cover

| Function                                               | What to mock                     | What to assert                                                             |
| ------------------------------------------------------ | -------------------------------- | -------------------------------------------------------------------------- |
| `api_utils.uniprot_accessions_to_uniparc_accessions()` | UniProt REST endpoint            | Correct `UID → UPI` mapping returned; `None` for UIDs absent from response |
| `structure_utils.download_alphafold_structure()`       | AlphaFold API + CIF download URL | CIF file written to expected path; returns `False` on `404`                |
| `structure_utils.fetch_structures()`                   | Both API layers                  | Files written with expected names for each UID in the input FASTA          |

---

## 6. Property-Based Tests

**Goal:** Generate random but valid inputs and assert that invariants always hold, catching edge cases not covered by handwritten examples.

**Location:** Alongside the unit tests they augment (e.g. `tests/unit/utils/test_fasta_utils.py`)

**Tooling:** `hypothesis`

```bash
pip install hypothesis
```

### Good candidates

| Function / class              | Invariant to assert                                                                      |
| ----------------------------- | ---------------------------------------------------------------------------------------- |
| `extract_uid_from_recordID()` | Any string of the form `sp\|XXXXX\|...` always returns the second pipe-delimited segment |
| `DataNormalizer` min-max      | Output is always in `[0, 1]` for any non-constant numeric input                          |
| Score DataFrame merging       | Merged output always has `≤ min(len(df1), len(df2))` rows for an inner join              |
| `ConfigFile` round-trip       | Any dict of string key-value pairs survives a write → read round-trip unchanged          |

---

## Recommended Directory Layout

```
tests/
├── __init__.py
├── conftest.py                       # Shared fixtures + e2e skip guards
├── fixtures/
│   ├── small.fasta
│   ├── small.clans
│   ├── small.tsv
│   ├── small.cif
│   └── small_roi.cif
├── unit/                             # Already implemented — see TESTING_GUIDE.md
│   ├── utils/
│   ├── core/
│   ├── similarity/
│   └── evaluation/                   # TODO: data_normalizer, cluster_analyzer
├── integration/
│   ├── test_clans_round_trip.py
│   ├── test_pipeline_mocked.py
│   └── test_scores_evaluator.py
├── e2e/
│   ├── test_pipeline_foldseek.py
│   └── test_pipeline_usalign.py
├── regression/
│   └── test_output_stability.py
└── contract/
    ├── test_api_utils_http.py
    └── test_structure_utils_http.py
```

---

## Priority Order

1. **Missing unit tests** (`evaluation/data_normalizer`, `evaluation/cluster_analyzer`) — self-contained, no I/O; complete the unit layer first
2. **Integration: CLANS round-trip** — uses real fixture files; highest value for catching regressions in the core data format
3. **Integration: pipeline with mocked tool** — validates the full composition of components without external dependencies
4. **Regression tests** — run after the integration layer is green; lock in stability before any larger refactor
5. **Contract tests** — decouple the test suite from live network access; unblock CI in offline environments
6. **E2E tests** — require real tool installations; run manually or in a dedicated CI environment
7. **Property-based tests** — add incrementally to the highest-risk pure functions

---

## Tooling Summary

| Library         | Purpose                       | Install                  |
| --------------- | ----------------------------- | ------------------------ |
| `pytest`        | Test runner                   | `pip install pytest`     |
| `pytest-cov`    | Coverage reporting            | `pip install pytest-cov` |
| `unittest.mock` | Mocking (stdlib)              | built-in                 |
| `responses`     | HTTP mocking                  | `pip install responses`  |
| `syrupy`        | Snapshot / regression testing | `pip install syrupy`     |
| `hypothesis`    | Property-based testing        | `pip install hypothesis` |

Run all non-E2E tests (unit + integration + regression + contract):

```bash
pytest tests/ -m "not e2e" --cov=src/clans3d --cov-report=term-missing
```
