# Issue #47 Plan: Optional User-Provided Structure Database (CIF)

## Scope

This plan addresses issue #47:

- Add support for using a user-provided database of local CIF files.
- When local structures are provided, skip AlphaFold downloads.
- Preserve existing behavior as default when no local structure database is provided.

Issue statement summary:

- Current pipeline always downloads structures from AlphaFold.
- Users should be able to provide their own CIF structure collection and bypass download.

## Goals

1. Add an opt-in path for local structure sourcing.
2. Keep backward compatibility for all current CLI/config usage.
3. Preserve data contracts used by downstream scoring and CLANS generation.
4. Ensure robust validation and clear error messaging for local structure mode.
5. Cover new behavior with unit and integration tests.

## Non-Goals

1. Supporting non-CIF structure formats (PDB, MMTF) in this issue.
2. Changing the similarity scoring contract or CLANS file format.
3. Implementing remote custom databases (HTTP/S3) in this issue.

## Current Architecture (Relevant)

1. CLI builds a PipelineConfig object.
2. Pipeline Step 1 calls fetch_structures(...) and always downloads from AlphaFold based on input IDs.
3. Returned contract from Step 1 is:
   - dict[str, tuple[int, int] | None]
   - uid to optional region
4. Downstream steps depend on structure files being present in structures_dir as:
   - <uid>.cif

Implication:

- The local database feature should preserve the same internal contract and filenames to avoid broad downstream changes.

## Proposed Design

### 1. New user-facing option

Add one optional argument:

- --structures_db <PATH>

Meaning:

- Path to a directory containing CIF files named <UniProtAccession>.cif.
- If set, pipeline uses local files instead of downloading from AlphaFold.

Behavior:

- Not set: keep current AlphaFold download behavior.
- Set: local structure mode is activated, download step is skipped.

### 2. Config model updates

Add a new PipelineConfig field:

- structures_db: str | None = None

Validation rules:

1. If structures_db is set, path must exist and be a directory.
2. Pipeline input validation still applies (input file exists, extension/type compatibility, tool checks).
3. Region constraints continue to apply exactly as today.

### 3. Structure preparation strategy

Introduce a new utility path in structure_utils:

- prepare_structures_from_local_db(input_file_path, input_file_type, structures_db, output_dir)

Responsibilities:

1. Parse input UIDs and optional regions using existing logic.
2. For each requested UID:
   - Expect structures_db/<uid>.cif.
   - If missing, warn and skip that UID (same spirit as download failure behavior).
3. Copy CIF files into output_dir as <uid>.cif for working isolation.
4. If region is provided, extract region into the copied output file (same as AlphaFold path).
5. Return dict[uid] = region for successfully prepared structures only.

Why copy into output_dir:

- Maintains existing pipeline/tool assumptions.
- Avoids mutating user-provided database files.
- Keeps run-to-run isolation and reproducibility.

### 4. Pipeline branching

In ClansPipeline.fetch_structures():

1. If config.structures_db is set:
   - Log that local database mode is active.
   - Call prepare_structures_from_local_db(...).
2. Else:
   - Keep existing fetch_structures(...) download path.
3. Keep existing RuntimeError if no structures are available after filtering.

### 5. CLI and config-file parity

Update CLI parser:

- Add --structures_db argument to main parser.

Ensure config file compatibility:

- Support -structures_db <PATH> in config files, automatically mapped to --structures_db by existing config loader.

### 6. Documentation updates

Update at minimum:

1. README.md
   - CLI options table.
   - Configuration file example.
   - Usage examples for local DB mode.
2. AGENTS.md
   - CLI argument table and examples.
   - Feature implementation checklist references.
3. docs/TESTING_GUIDE.md
   - If new test files/classes are added.

## Detailed Implementation Plan (Step by Step)

### Step 0: Baseline and contract lock

1. Confirm current fetch_structures and pipeline assumptions with tests.
2. Record key invariants:
   - structures_dir contains <uid>.cif
   - returned mapping only includes successfully prepared structures

Exit criteria:

- Invariants are documented in tests/comments where needed.

### Step 1: Add CLI + config plumbing

1. Add --structures_db to core/cli.py.
2. Pass args.structures_db into PipelineConfig in main.py.
3. Add structures_db field and validation in core/pipeline.py.

Exit criteria:

- Parsing works from CLI and config file.
- Invalid path produces clear ValueError.

### Step 2: Implement local-DB structure preparation utility

1. Add helper(s) in utils/structure_utils.py to:
   - parse requested UID/region tasks from input
   - copy matching CIF files from structures_db
   - apply optional region extraction on copied files
2. Reuse existing region extraction function for consistency.
3. Keep warning-level logs for missing/invalid entries and continue processing.

Exit criteria:

- Utility returns expected uid-to-region mapping.
- Output directory contains only successfully prepared structures.

### Step 3: Integrate pipeline branching

1. Update ClansPipeline.fetch_structures() to branch on structures_db.
2. Preserve old behavior when structures_db is not provided.
3. Preserve failure behavior when zero structures are available.

Exit criteria:

- Existing workflows unchanged.
- Local DB mode bypasses AlphaFold network path completely.

### Step 4: Testing

Unit tests:

1. tests/unit/core/test_cli.py
   - parse --structures_db
   - config override behavior remains correct
2. tests/unit/core/test_pipeline_validate.py
   - valid directory accepted
   - missing/non-directory path rejected
3. tests/unit/utils/test_structure_utils.py
   - local DB utility returns only successful UIDs
   - missing CIF logs warning and skips UID
   - region extraction is applied to copied file
   - no mutation of source structures_db files

Integration tests:

1. tests/integration/test_pipeline_mocked.py (or new file)
   - Local DB mode path with no network calls.
   - Verify pipeline.run() output with local structures.
   - Assert AlphaFold download function is not called in local mode.

Regression checks:

1. Existing FASTA/TSV behavior without structures_db remains green.

Exit criteria:

- New tests pass.
- Existing relevant test suites pass.

### Step 5: Documentation and examples

1. Update README examples with local DB usage.
2. Update AGENTS CLI table and examples.
3. Add a small example config snippet showing -structures_db.

Exit criteria:

- User-facing docs fully describe both modes.

## Validation and Error-Handling Rules

1. structures_db path invalid:
   - raise ValueError with actionable message.
2. UID requested but CIF not found in local DB:
   - log warning and skip UID.
3. Region extraction fails for one structure:
   - log warning and skip UID (do not abort whole run).
4. No structures available after filtering:
   - raise RuntimeError, consistent with current behavior.

## Risks and Mitigations

### Risk 1: UID-to-filename mismatches in user DB

Mitigation:

- Explicitly document required naming convention <uid>.cif.
- Log missing file warnings with UID and expected path.

### Risk 2: Local DB mode accidentally mutates source files

Mitigation:

- Always copy source CIF to working structures_dir first.
- Perform region extraction only on copied file.

### Risk 3: Regressions in existing AlphaFold flow

Mitigation:

- Keep branch minimal and isolated to structure preparation step.
- Preserve existing code path untouched when structures_db is None.

### Risk 4: Silent partial dataset due to many missing CIFs

Mitigation:

- Keep aggregated success/failure logging counts.
- Maintain hard failure when successful set is empty.

## Definition of Done

1. User can run pipeline with --structures_db and no AlphaFold downloads occur.
2. Existing workflow without --structures_db behaves exactly as before.
3. Local DB path supports FASTA/A2M/A3M/TSV input UID parsing and optional regions.
4. Unit and integration tests cover success, skip, and failure cases.
5. README and AGENTS docs are updated with new option and examples.

## Suggested Implementation Order

1. CLI/config plumbing and PipelineConfig validation.
2. Local DB utility implementation in structure_utils.
3. Pipeline branching integration.
4. Test additions and updates.
5. Documentation updates.

This order minimizes regression risk by introducing typed/configured entry points first, then adding behavior in isolated utilities, then wiring orchestration, then validating and documenting.
