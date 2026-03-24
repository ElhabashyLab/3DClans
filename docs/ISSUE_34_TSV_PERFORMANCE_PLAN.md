# Issue #34 Plan: TSV Path Sequential Network Bottlenecks

## Scope

This plan targets issue #34:

- TSV workflow has sequential API bottlenecks in UniProt to UniParc mapping and FASTA retrieval.
- Goal is to reduce TSV wall-clock time while preserving output semantics, deterministic ordering, and error handling behavior.

## Status Update (2026-03-22)

### Latest benchmark snapshot (1500 TSV input)

Source: `benchmark_output/run_20260322_174112/benchmark_results.csv`

| Tool     | Score Type | Num Structures | Download Time (s) | Computation Time (s) | Generation Time (s) | Total Time (s) | Num Scores |
| -------- | ---------- | -------------- | ----------------- | -------------------- | ------------------- | -------------- | ---------- |
| foldseek | evalue     | 1499           | 1299.71           | 44.75                | 33.60               | 1378.05        | 837295     |
| foldseek | TM         | 1499           | 1299.71           | 76.98                | 33.30               | 1409.99        | 837295     |
| USalign  | -          | 1499           | 1299.71           | 4082.37              | 44.88               | 5426.95        | 1122751    |

Interpretation notes:

- In the benchmark implementation, `Download Time` includes both structure fetching and cleaned FASTA preparation during benchmark initialization.
- For foldseek runs, this initialization stage is currently the dominant cost.
- For USalign runs, computation remains the dominant cost.

### What is done

- Step 1 implemented: `generate_fasta_from_uids_with_regions` now performs bounded concurrent FASTA retrieval while preserving deterministic output ordering.
- Step 3 implemented as decided: existing `-w/--workers` (`PipelineConfig.download_workers`) is reused for TSV FASTA retrieval; default remains 10.
- Step 5 partially implemented:
  - Added/updated unit tests for deterministic ordering under concurrency.
  - Added worker validation tests.
  - Added pipeline test coverage for forwarding `download_workers` into TSV FASTA generation.

### What is left

- Step 0 (remaining): add stage-specific timing around mapping and FASTA retrieval paths to isolate TSV API costs from structure-download costs in benchmarks.
- Step 2: optimize UniProt mapping orchestration (batch overlap, bounded concurrency, preserve mapping contract and determinism).
- Step 4: improve observability for long TSV runs (explicit counters/progress for mapping batches and FASTA retrieval).
- Step 5 (remaining): add targeted integration/regression coverage for Step 2 behavior and benchmark comparison artifacts.

## 1. Problem Breakdown

### Problem A: Blocking mapping lifecycle per batch

Current flow in `uniprot_accessions_to_uniparc_accessions` does this per batch:

1. Submit mapping job.
2. Poll until done.
3. Fetch all pages.
4. Move to next batch.

Why this hurts:

- Submit/poll/fetch is serialized across batches.
- Network idle time accumulates while waiting for each batch completion.

### Problem B: Sequential pagination fetch

`_fetch_paginated_results` follows next links in a single-thread loop.

Why this hurts:

- Even when each page is small, network latency stacks linearly.
- End-to-end mapping stage grows with number of pages.

Note:

- This step is naturally order-dependent because next URL is discovered from the previous response.
- Biggest practical gains come from reducing timeouts/retries overhead and overlapping at higher levels, not from unsafe parallel page fetches.

### Problem C: Per-UID sequential FASTA retrieval

`generate_fasta_from_uids_with_regions` downloads one sequence at a time by calling `download_fasta_record` in a loop.

Why this hurts:

- End-to-end time is roughly sum of per-request latencies.
- Large TSV inputs are dominated by request round-trip accumulation.

### Problem D: No concurrency controls specific to TSV FASTA stage

There is no worker setting for this stage.

Why this hurts:

- Cannot tune throughput based on network/API constraints.
- Risk of accidental overloading if concurrency is added without limits.

### Problem E: Output ordering and fallback semantics must not regress

Current behavior guarantees a record per input UID order, using mock-up records on failures.

Why this matters:

- Parallelization can reorder completion times.
- Reordering outputs or changing fallback behavior would be a regression.

## 2. Solution Strategy ( Step by Step)

### Step 0: Establish baseline timing and correctness

1. Add temporary timing logs around:
   - mapping conversion stage,
   - FASTA download stage,
   - total `generate_fasta_from_uids_with_regions` time.
2. Run on at least one medium and one large TSV file.
3. Save baseline metrics in issue notes.

Exit criteria:

- Baseline stage timings and output hashes/record counts are captured.

### Step 1: Add bounded concurrency for FASTA retrieval (highest impact in issue #34)

1. Introduce worker-based concurrent retrieval in `generate_fasta_from_uids_with_regions`.
2. Use a thread pool with a conservative default (for example 8 to 16 workers).
3. Keep deterministic output order by pre-allocating results by input index and writing in original UID order.
4. Preserve current fallback behavior:
   - on ValueError or RequestException, create mock-up record.

Implementation notes:

- Keep `download_fasta_record` as single-record primitive.
- Add a small internal worker helper for one UID task.

Exit criteria:

- Same output order and failure semantics as before.
- Clear reduction in TSV FASTA stage time.

### Step 2: Improve mapping batch orchestration without changing data contract

1. Keep page traversal sequential per job (safe and simple).
2. Overlap at batch level:
   - submit multiple jobs first (bounded),
   - then poll/fetch completed jobs with bounded worker pool.
3. Merge per-batch mappings into final mapping while filling unmapped entries with None.

Implementation notes:

- Preserve `dict[str, str | None]` contract.
- Keep deterministic filling by iterating original accessions when constructing the final map.

Exit criteria:

- Mapping stage latency decreases on multi-batch inputs.
- Unmapped behavior remains unchanged.

### Step 3: Reuse existing workers flag for TSV FASTA stage

1. Reuse existing `-w/--workers` (`PipelineConfig.download_workers`) for TSV FASTA concurrency.
2. Keep default at 10 workers.
3. Do not add new CLI/config options for this stage unless profiling shows a need.

Exit criteria:

- Existing workers setting controls throughput consistently across structure download and TSV FASTA retrieval.
- Default of 10 remains stable for typical runs.

### Step 4: Harden long-run observability

1. Add progress logs for FASTA retrieval:
   - processed/total,
   - success/failure counters.
2. Add progress logs for mapping batches:
   - submitted/completed batches,
   - per-stage elapsed time.

Exit criteria:

- Long TSV runs provide actionable progress visibility.

### Step 5: Validation and regression safety

1. Unit tests:
   - keep output order deterministic under concurrency,
   - verify mock-up fallback per failed UID,
   - verify mapping contract with missing IDs.
2. Integration tests:
   - TSV workflow with patched network/subprocess,
   - compare generated FASTA contents before/after change.
3. Performance validation:
   - measure stage timings from Step 0 and compare.

Exit criteria:

- No semantic regressions.
- Measurable runtime improvement for TSV path.

## 3. Proposed Implementation Order

1. Step 1 (FASTA concurrency).
2. Step 5 partial (tests for Step 1).
3. Step 2 (mapping orchestration).
4. Step 5 partial (tests for Step 2).
5. Step 3 (tuning knobs).
6. Step 4 (logging polish).
7. Final benchmark pass and documentation updates.

Reasoning:

- Step 1 is direct and high-yield.
- Step 2 is second-largest issue #34 lever.
- Config and logging come after correctness/performance gains are proven.

## 4. Risks and Mitigations

### Risk: API rate limits or transient errors increase under concurrency

Mitigation:

- Bounded workers with conservative defaults.
- Reuse issue #35 follow-up for shared session/retry policy.

### Risk: Nondeterministic output order

Mitigation:

- Store results by original input index and serialize in input order.

### Risk: Hidden behavior change for missing IDs

Mitigation:

- Snapshot tests of generated FASTA IDs and sequence placeholders.

## 5. Definition of Done for Issue #34

- TSV end-to-end runtime improves on medium and large datasets.
- No change in output ordering or fallback semantics.
- Mapping and FASTA stages show reduced wall-clock time.
- Unit/integration tests cover concurrency and failure paths.
- Logging clearly reports progress for long TSV runs.
