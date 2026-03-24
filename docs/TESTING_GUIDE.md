# Testing Guide for 3DClans

All commands below assume you are in the project root and the virtual environment is active.

---

## Running Tests

### Run all tests (unit + integration + regression + e2e when selected)

```bash
pytest tests/
```

### Run all tests except e2e

```bash
pytest tests/ -m "not e2e"
```

### Run all unit tests

```bash
pytest tests/unit/
```

### Run with coverage report

```bash
pytest tests/unit/ --cov=src/clans3d --cov-report=term-missing
```

### Target a specific module

```bash
pytest tests/unit/utils/test_fasta_utils.py
```

### Run a single test by name

```bash
pytest tests/unit/utils/test_fasta_utils.py::TestExtractUidFromRecordID
```

---

## Unit Test Coverage

Unit tests live in `tests/unit/` and are organized to mirror the `src/clans3d/` package layout. All external dependencies (subprocesses, network, filesystem) are mocked.

### `utils/fasta_utils.py` — `tests/unit/utils/test_fasta_utils.py`

| Covered class / method                   | What is tested                                                                                         |
| ---------------------------------------- | ------------------------------------------------------------------------------------------------------ |
| `TestExtractUidFromRecordID`             | Standard `sp\|UID\|...`, `tr\|UID\|...`, plain UID, version suffix, region suffix, underscore-style    |
| `TestExtractRegionFromRecord`            | Header with and without region, large region numbers                                                   |
| `TestAddRegionToRecord`                  | Valid region trims sequence and updates ID; `None` leaves unchanged; out-of-bounds raises `ValueError` |
| `TestCreateMockUpRecord`                 | Correct ID format with/without region; sequence is `"not_found"`                                       |
| `TestExtractUidsFromFasta`               | Reads a fixture FASTA and returns UIDs in order                                                        |
| `TestExtractRecordsFromFasta`            | Returns correct number of `SeqRecord` objects                                                          |
| `TestCleanFastaFile`                     | Retains only records in the allow-list; writes output to file                                          |
| `TestCopyRecordsFromFasta`               | Only records matching given UIDs are written                                                           |
| `TestDownloadFastaRecord`                | Correct record is fetched from a mocked HTTP endpoint                                                  |
| `TestRemoveNonExistingUniprotAccessions` | Filters out accessions absent from a mocked UniProt response                                           |
| `TestGenerateFastaFromUidsWithRegions`   | Generates FASTA lines with region annotations from a UID/region mapping                                |
| `TestCleanAlignedSequence`               | Removes alignment gap chars (`.` and `-`), preserves lowercase insertions, uppercases output           |
| `TestGenerateFastaFromAlignmentFile`     | Generates cleaned FASTA from A2M/A3M-style alignments, handles invalid headers and empty inputs         |
| `TestCleanedSequenceMatchesRegionLength` | Ensures cleaned sequence lengths match requested region spans across representative cases                |

### `utils/api_utils.py` — `tests/unit/utils/test_api_utils.py`

| Covered class / method      | What is tested                                                           |
| --------------------------- | ------------------------------------------------------------------------ |
| `TestChunked`               | Even split, uneven split, empty list, chunk larger than list             |
| `TestSubmitIdmappingJob`    | Mocked HTTP POST; correct job ID returned; error response raises         |
| `TestWaitForResults`        | Polls until status is `FINISHED`; raises on timeout or server error      |
| `TestFetchPaginatedResults` | Iterates pages; assembles full result list; handles single-page response |

### `utils/dependency_checks.py` — `tests/unit/utils/test_dependency_checks.py`

| Covered class / method       | What is tested                                                                 |
| ---------------------------- | ------------------------------------------------------------------------------ |
| `TestCheckExternalTool`      | Returns path when `shutil.which` finds the tool; returns `None` when not found |
| `TestVerifyToolDependencies` | Raises `RuntimeError` when tool missing; does not raise when present           |

### `utils/file_utils.py` — `tests/unit/utils/test_file_utils.py`

| Covered class / method | What is tested                                                                      |
| ---------------------- | ----------------------------------------------------------------------------------- |
| `TestResetDirContent`  | Creates directory when missing; removes existing files before recreation            |
| `TestCopyDirContent`   | Files and nested subdirectories are copied; target created when missing             |
| `TestDownloadFile`     | Returns `True` and writes file on HTTP 200; returns `False` and warns on HTTP error |

### `utils/log.py` — `tests/unit/utils/test_log.py`

| Covered class / method | What is tested                                                                                                    |
| ---------------------- | ----------------------------------------------------------------------------------------------------------------- |
| `TestSetupLogging`     | `verbose=True` → `DEBUG`; `quiet=True` → `ERROR`; default → `INFO`; calling twice does not add duplicate handlers |

### `utils/structure_utils.py` — `tests/unit/utils/test_structure_utils.py`

| Covered class / method           | What is tested                                                                                            |
| -------------------------------- | --------------------------------------------------------------------------------------------------------- |
| `TestLogInterval`                | Correct log interval for small, medium, and large totals                                                  |
| `TestParseRegionFromTsvRow`      | Valid TSV row returns `(start, end)`; malformed row raises                                                |
| `TestFetchAlphafoldCifUrl`       | Mocked API returns correct CIF URL; missing entry returns `None`                                          |
| `TestDownloadAlphafoldStructure` | File written on success; returns `False` on API error; region extraction called only when region provided |
| `TestFetchStructures`            | Dispatches to correct handler for FASTA/A2M/TSV; raises `ValueError` for unsupported type                 |
| `TestRunDownloadTasks`           | Parallel download task behavior: successful-result filtering, worker count usage, exception propagation   |
| `TestExtractRegionOfProtein`     | Output CIF contains only residues within the specified range                                              |

### `core/clans_file.py` — `tests/unit/core/test_clans_file.py`

| Covered class / method      | What is tested                                                                                  |
| --------------------------- | ----------------------------------------------------------------------------------------------- |
| `TestClansFileConstruction` | Raises `ValueError` when both `fasta_records` and `path_to_fasta` are `None`; accepts either    |
| `TestValidation`            | Raises when coordinate count ≠ sequence count; raises when coordinate/score UIDs not in records |
| `TestClansFileStr`          | Output contains `<param>`, `<seq>`, `<pos>`, `<hsp>` sections; indices match FASTA record order |

### `core/clans_file_generator.py` — `tests/unit/core/test_clans_file_generator.py`

| Covered class / method          | What is tested                                                                 |
| ------------------------------- | ------------------------------------------------------------------------------ |
| `TestExtractBlock`              | Returns lines between tags; returns `None` when tag missing                    |
| `TestParseNumberOfSequences`    | Correct count parsed; raises `ValueError` when `sequences=` line is absent     |
| `TestParseParamBlock`           | Returns dict of params; handles absent param block                             |
| `TestParseParamBlockEdgeCases`  | Empty block, duplicate keys, whitespace handling                               |
| `TestParseFastaBlock`           | Returns correct number of `SeqRecord` objects; raises on missing `<seq>` block |
| `TestParsePosBlock`             | Returns correct `(idx, x, y, z)` tuples; raises on malformed lines             |
| `TestParseScoresBlock`          | Returns DataFrame with correct columns and dtypes                              |
| `TestParseScoresBlockEmpty`     | Returns empty DataFrame for empty `<hsp>` block                                |
| `TestNormalizeScoresFormat`     | Removes duplicate A↔B/B↔A pairs; result is smaller than input                  |
| `TestGenerateRandomCoordinates` | Returns one tuple per UID; x/y/z all in `[0, 1]`                               |
| `TestGenerateClansFile`         | Generated file is written and contains correct sequence count                  |

### `core/config_file.py` — `tests/unit/core/test_config_file.py`

| Covered class / method   | What is tested                                                                                  |
| ------------------------ | ----------------------------------------------------------------------------------------------- |
| `TestWriteAndReadConfig` | Written file contains `-key value` lines; comments ignored on read; round-trip preserves values |
| `TestConfigToArgv`       | Returns `["--key", "value", ...]` list matching config contents                                 |

### `core/enums` — `tests/unit/core/test_enums.py`

| Covered class / method | What is tested                                                                     |
| ---------------------- | ---------------------------------------------------------------------------------- |
| `TestInputFileType`    | `InputFileType("fasta") == InputFileType.FASTA`; invalid value raises `ValueError` |
| `TestToolType`         | `ToolType("foldseek") == ToolType.FOLDSEEK`; invalid value raises `ValueError`     |

### `core/pipeline.py` — `tests/unit/core/test_pipeline_validate.py`

| Covered class / method       | What is tested                                                                                                       |
| ---------------------------- | -------------------------------------------------------------------------------------------------------------------- |
| `TestPipelineConfigDefaults` | `structures_dir`, `output_dir`, `cleaned_input_storage` use expected default paths                                   |
| `TestClansPipelineValidate`  | Raises `FileNotFoundError` for missing input; raises `ValueError` on extension mismatch and invalid tool/score combo |
| `TestFetchStructures`        | Correct structures-fetch function called for the configured input type                                               |

### `core/cli.py` — `tests/unit/core/test_cli.py`

| Covered class / method | What is tested                                                                                               |
| ---------------------- | ------------------------------------------------------------------------------------------------------------ |
| `TestParseArgs`        | Required args → correct namespace; config file merged and overridden by CLI; invalid tool/type → parse error |

### `similarity/foldseek.py` — `tests/unit/similarity/test_foldseek.py`

| Covered class / method  | What is tested                                                                             |
| ----------------------- | ------------------------------------------------------------------------------------------ |
| `TestGetOutputColumns`  | Returns `"evalue"` columns for `"evalue"`; TM columns for `"TM"`; raises for unknown score |
| `TestCleanScoresEvalue` | Self-hits removed; symmetric duplicates removed; `score` column equals `evalue`            |
| `TestCleanScoresTM`     | `score = 1 - max(TM1, TM2)`; output has `score` column, no `TM1`/`TM2`                     |
| `TestDetectPhase`       | Returns correct label when known keyword in chunks; returns `""` when no match             |
| `TestCreateDatabase`    | Subprocess called with expected arguments; working directory set correctly                 |
| `TestParseOutputEvalue` | Mocked output produces DataFrame with correct evalue-based scores                          |
| `TestParseOutputTM`     | Mocked output produces DataFrame with correct TM-based scores                              |

### `similarity/usalign.py` — `tests/unit/similarity/test_usalign.py`

| Covered class / method | What is tested                                                                                                |
| ---------------------- | ------------------------------------------------------------------------------------------------------------- |
| `TestCountDataLines`   | Ignores `#` and `Warning` lines; counts remaining non-empty lines                                             |
| `TestParseOutput`      | Mocked TSV output → DataFrame with correct columns; `score = 1 - max(TM1, TM2)`; extensions stripped from IDs |

### `similarity/struct_sim_computer.py` — `tests/unit/similarity/test_struct_sim_computer.py`

| Covered class / method | What is tested                                                                                        |
| ---------------------- | ----------------------------------------------------------------------------------------------------- |
| `TestCreateTool`       | Returns `Foldseek` for `ToolType.FOLDSEEK`; `USalign` for `ToolType.USALIGN`; raises for unknown type |
| `TestRun`              | Logs warning when returned score count ≠ expected                                                     |

### `similarity/struct_sim_tool.py` — `tests/unit/similarity/test_struct_sim_tool.py`

| Covered class / method | What is tested                                                   |
| ---------------------- | ---------------------------------------------------------------- |
| `TestExecuteRun`       | Subprocess called with correct arguments; stdout/stderr captured |

---

## Tooling Summary

| Library         | Purpose            | Install                  |
| --------------- | ------------------ | ------------------------ |
| `pytest`        | Test runner        | `pip install pytest`     |
| `pytest-cov`    | Coverage reporting | `pip install pytest-cov` |
| `unittest.mock` | Mocking (stdlib)   | built-in                 |

---

## Integration Tests

Integration tests live in `tests/integration/` and verify that pipeline components work together end-to-end using real fixture files. External dependencies (AlphaFold downloads, Foldseek/USalign subprocesses) are patched out.

### Run

```bash
pytest tests/integration/
```

### `tests/integration/test_a2m_sequence_cleaning.py`

Validates region-based cleaning against known UniProt regions and confirms lowercase residues are treated as real amino acids (not removed as gaps).

| Class                               | What is verified                                                                 |
| ----------------------------------- | -------------------------------------------------------------------------------- |
| `TestCleanedSequenceMatchesUniProt` | Cleaned extracted sequence matches expected UniProt region content and length    |
| `TestRealA2MFileAgainstUniProt`     | Real fixture A2M entries match expected UniProt-derived cleaned sequences         |
| `TestLowercaseAreRealResidues`      | Lowercase letters are preserved as residues; removing them would lose information |

### `tests/integration/test_pipeline_mocked.py` — `TestPipelineMocked` (FASTA input)

Uses `tests/fixtures/small.fasta` (3 sequences) with `fetch_structures` and `StructSimComputer` mocked to return pre-built data. Runs `ClansPipeline.run()` end-to-end and asserts:

| Test                                             | What is verified                                        |
| ------------------------------------------------ | ------------------------------------------------------- |
| `test_run_writes_clans_file`                     | Output `.clans` file is created on disk                 |
| `test_run_output_has_correct_sequence_count`     | `sequences=3` appears in the output                     |
| `test_run_output_contains_all_hsp_entries`       | All 3 pairwise score entries appear in `<hsp>`          |
| `test_run_returns_existing_fasta_path`           | Cleaned FASTA path is returned and the file exists      |
| `test_run_output_contains_all_required_sections` | All required CLANS sections are present                 |
| `test_run_raises_when_no_structures_available`   | `RuntimeError` raised when no structures are downloaded |

### `tests/integration/test_pipeline_mocked.py` — `TestPipelineMockedTSV` (TSV input)

Uses `tests/fixtures/small.tsv` with the same mocks as above, plus `uniprot_accessions_to_uniparc_accessions` and `download_fasta_record` patched to avoid network calls. Exercises the TSV code path where cleaned FASTA is built from scratch rather than copied from the input file.

| Test                                         | What is verified                                   |
| -------------------------------------------- | -------------------------------------------------- |
| `test_run_writes_clans_file`                 | Output `.clans` file is created on disk            |
| `test_run_output_has_correct_sequence_count` | `sequences=3` appears in the output                |
| `test_run_output_contains_all_hsp_entries`   | All 3 pairwise score entries appear in `<hsp>`     |
| `test_run_returns_existing_fasta_path`       | Cleaned FASTA path is returned and the file exists |

---

## Regression Tests

Regression tests live in `tests/regression/` and lock in the structural properties of generated CLANS files so that refactors cannot silently break output shape.

### Run

```bash
pytest tests/regression/
```

### `tests/regression/test_clans_output.py` — `TestClansOutputStructure`

Calls `ClansFileGenerator.generate_clans_file()` directly with `tests/fixtures/small.fasta` and fixed mock scores. No network or subprocess calls.

| Test                                 | What is locked in                                |
| ------------------------------------ | ------------------------------------------------ |
| `test_sequence_count_stable`         | `sequences=3` in output                          |
| `test_coordinate_count_stable`       | Exactly 3 coordinate lines in `<pos>`            |
| `test_score_entry_count_stable`      | Exactly 3 score lines in `<hsp>`                 |
| `test_all_required_sections_present` | All section tags are present                     |
| `test_score_format_stable`           | Each score line matches `idx1 idx2:score` format |
| `test_fasta_sequences_preserved`     | All three UIDs appear in `<seq>`                 |

---

## End-to-End (E2E) Tests

E2E tests live in `tests/e2e/` and require Foldseek and/or USalign to be installed in `PATH`. They download real AlphaFold structures for the sequences in `examples/small_fasta_files/5.fasta` and run the full pipeline. Tests are skipped automatically when the required binary is absent and must be requested explicitly via the `e2e` mark.

### Run

```bash
pytest tests/e2e/ -m e2e
```

### `tests/e2e/test_pipeline_e2e.py`

FASTA input: `examples/small_fasta_files/5.fasta` (4 real MYOD1 orthologs + 1 unavailable entry; ≥ 3 structures expected to download successfully).
TSV input: `examples/small_tsv_files/5.tsv` (same accessions in TSV format; sequences fetched from UniProt).

| Class                        | Tool                | Input | Tests                                                              |
| ---------------------------- | ------------------- | ----- | ------------------------------------------------------------------ |
| `TestPipelineFoldseekEvalue` | Foldseek (evalue)   | FASTA | File created, ≥ 3 sequences, ≥ 1 score entry, all sections present |
| `TestPipelineFoldseekTM`     | Foldseek (TM score) | FASTA | File created, all scores in `[0, 1]`                               |
| `TestPipelineUSalign`        | USalign             | FASTA | File created, ≥ 3 sequences, ≥ 1 score entry                       |
| `TestCliEntrypoint`          | Foldseek            | FASTA | `3dclans` subprocess exits with code `0`                           |
| `TestPipelineTSVInput`       | Foldseek (evalue)   | TSV   | File created, ≥ 3 sequences, ≥ 1 score entry, cleaned FASTA written |
