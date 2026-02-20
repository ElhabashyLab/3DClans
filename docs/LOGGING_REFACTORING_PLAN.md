# Logging & Console Output Refactoring Plan

**Related issue:** [#30 — Provide information about progress](https://github.com/ElhabashyLab/Clans-3D/issues/30)  
**Scope:** All modules except `evaluation/` and `dataset_generator/`

---

## 1. Problem Statement

### 1.1 No unified output system

All console output uses raw `print()` calls. There is no logging framework, no consistent format, and no way to control verbosity. The 48 `print()` statements across the in-scope modules fall into several categories but are visually indistinguishable from each other.

### 1.2 No progress feedback (Issue #30)

Long-running operations (structure downloads, similarity computation, CLANS file generation) provide little or no progress information. The user has no visibility into how far along the process is.

### 1.3 Inconsistent message styles

- Some messages use `\n` prefixes, some don't
- Warnings are printed as plain text (`"Warning: ..."`) — not distinguishable from info
- Errors mix `print()` with exceptions
- No timestamps, no severity indicators, no module attribution

---

## 2. Current Print Statement Audit

### 2.1 `main.py` (2 prints)

| Line | Message                                    | Category   |
| ---- | ------------------------------------------ | ---------- |
| 22   | `CLANS file: {clans_file_path}`            | **Result** |
| 23   | `Cleaned FASTA: {cleaned_input_file_path}` | **Result** |

### 2.2 `core/clans_file_generator.py` (3 prints)

| Line | Message                                   | Category          |
| ---- | ----------------------------------------- | ----------------- |
| 26   | `Parsing CLANS file {clans_file_path}...` | **Step start**    |
| 184  | `Generating CLANS file {out_path}...`     | **Step start**    |
| 199  | `CLANS file generated at {out_path}`      | **Step complete** |

### 2.3 `similarity/struct_sim_computer.py` (3 prints)

| Line | Message                                                         | Category          |
| ---- | --------------------------------------------------------------- | ----------------- |
| 26   | `Computing structural similarity with {tool.name}...`           | **Step start**    |
| 35   | `Warning: {tool.name} did not return expected scores...`        | **Warning**       |
| 36   | `Structural similarity computation with {tool.name} completed.` | **Step complete** |

### 2.4 `similarity/struct_sim_tool.py` (1 print)

| Line | Message                                              | Category  |
| ---- | ---------------------------------------------------- | --------- |
| 37   | `Error running {self.name} with {self.command}: {e}` | **Error** |

### 2.5 `similarity/foldseek.py` (6 prints)

| Line | Message                                                     | Category  |
| ---- | ----------------------------------------------------------- | --------- |
| 80   | `Error running {self.name} with {create_db_command}: {e}`   | **Error** |
| 81   | `stdout: {e.stdout}`                                        | **Debug** |
| 82   | `stderr: {e.stderr}`                                        | **Debug** |
| 103  | `Error running {self.name} with {convertalis_command}: {e}` | **Error** |
| 104  | `stdout: {e.stdout}`                                        | **Debug** |
| 105  | `stderr: {e.stderr}`                                        | **Debug** |

### 2.6 `similarity/tmalign.py` (2 prints)

| Line | Message                                       | Category  |
| ---- | --------------------------------------------- | --------- |
| 58   | `Failed to parse output: DataFrame is empty.` | **Error** |
| 61   | `Error parsing output: {e}`                   | **Error** |

### 2.7 `similarity/usalign.py` (2 prints)

| Line | Message                                       | Category  |
| ---- | --------------------------------------------- | --------- |
| 56   | `Failed to parse output: DataFrame is empty.` | **Error** |
| 59   | `Error parsing output: {e}`                   | **Error** |

### 2.8 `utils/structure_utils.py` (7 prints)

| Line | Message                                               | Category                    |
| ---- | ----------------------------------------------------- | --------------------------- |
| 29   | `Downloading structure files in "{output_dir}"...`    | **Step start**              |
| 61   | `Downloaded {n} from {total} PDB files successfully.` | **Step complete / Summary** |
| 93   | `Downloaded {n} from {total} PDB files successfully.` | **Step complete / Summary** |
| 128  | `Failed to download {id} from {api_url}`              | **Warning** (non-fatal)     |
| 131  | `Failed to download {id} from {api_url}`              | **Warning** (non-fatal)     |
| 140  | `Failed to download {id} from {api_url}`              | **Warning** (non-fatal)     |
| 146  | `Failed to download {id} from {api_url}`              | **Warning** (non-fatal)     |

### 2.9 `utils/file_utils.py` (1 print)

| Line | Message                              | Category                |
| ---- | ------------------------------------ | ----------------------- |
| 21   | `Failed to download {url}: {str(e)}` | **Warning** (non-fatal) |

### 2.10 `utils/fasta_utils.py` (3 prints)

| Line    | Message                                      | Category       |
| ------- | -------------------------------------------- | -------------- |
| 205     | `Could not retrieve sequence for UID: {uid}` | **Warning**    |
| 220     | `Filtering...`                               | **Step start** |
| 222/227 | `Extracting uids ...`                        | **Step start** |

### 2.11 `utils/dependency_checks.py` (1 print)

| Line | Message                             | Category           |
| ---- | ----------------------------------- | ------------------ |
| 45   | `Found {tool_name} at: {tool_path}` | **Info** (verbose) |

### 2.12 `legacy/utils_old_clans.py` (1 print)

| Line | Message                    | Category  |
| ---- | -------------------------- | --------- |
| 24   | `running with {arguments}` | **Debug** |

### 2.13 `benchmark/benchmark.py` (26 prints)

Benchmark has its own formatted output with `=`/`-` delimiters and structured result tables. This is a dedicated reporting module and should retain its own style, but should migrate to the same logging infrastructure.

### Summary by category

| Category                         | Count | Shown at default verbosity? |
| -------------------------------- | ----- | --------------------------- |
| **Result** (final output paths)  | 2     | Always                      |
| **Step start/complete**          | 9     | Always                      |
| **Warning** (non-fatal failures) | 8     | Always                      |
| **Error** (fatal/tool failures)  | 7     | Always                      |
| **Debug** (stdout/stderr dumps)  | 4     | Verbose only                |
| **Info** (nice-to-know)          | 2     | Verbose only                |
| **Benchmark reporting**          | 26    | Own context                 |

---

## 3. Proposed Solution

### 3.1 Introduce Python `logging` module

Create a central logger configuration in a new module `src/clans3d/utils/log.py`:

```python
# src/clans3d/utils/log.py
import logging

def setup_logging(verbose: bool = False) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    formatter = logging.Formatter(
        fmt="[%(levelname)-7s] %(message)s"
    )
    handler = logging.StreamHandler()
    handler.setFormatter(formatter)

    root = logging.getLogger("clans3d")
    root.setLevel(level)
    root.addHandler(handler)
    root.propagate = False
```

Each module gets its own logger:

```python
import logging
logger = logging.getLogger(__name__)
```

### 3.2 Output format specification

**Default mode** (`-v` not set, level `INFO`):

```
[INFO   ] Found foldseek at: /usr/local/bin/foldseek
[INFO   ] Downloading structures... (5 proteins)
[WARNING] Failed to download P12345 — skipping
[INFO   ] Downloaded 4/5 PDB files
[INFO   ] Computing structural similarity with foldseek...
[INFO   ] Similarity computation completed (4 scores)
[INFO   ] Generating CLANS file output/5_cleaned.clans...
[INFO   ] CLANS file generated at output/5_cleaned.clans

Results:
  CLANS file:    output/clans_files/5_cleaned.clans
  Cleaned FASTA: work/input_file_storage/5_cleaned.fasta
```

**Verbose mode** (`-v` set, level `DEBUG`):
Adds debug-level details like subprocess commands, stdout/stderr from tools, individual download attempts, etc.

### 3.3 Level mapping for existing prints

| Current pattern                                         | New level                         | Visible by default? |
| ------------------------------------------------------- | --------------------------------- | ------------------- |
| Final result paths                                      | `INFO` (or plain print, see §3.5) | Yes                 |
| `"Downloading..."`, `"Computing..."`, `"Generating..."` | `INFO`                            | Yes                 |
| `"Downloaded X/Y..."`, `"completed."`                   | `INFO`                            | Yes                 |
| `"Warning: ..."`                                        | `WARNING`                         | Yes                 |
| `"Failed to download {id}"` (non-fatal, skipped)        | `WARNING`                         | Yes                 |
| `"Error running ..."` (tool failure)                    | `ERROR`                           | Yes                 |
| `stdout: ...`, `stderr: ...`                            | `DEBUG`                           | No — verbose only   |
| `"Found {tool} at: {path}"`                             | `DEBUG`                           | No — verbose only   |
| `"running with {args}"`                                 | `DEBUG`                           | No — verbose only   |
| `"Filtering..."`, `"Extracting uids..."`                | `DEBUG`                           | No — verbose only   |

### 3.4 Add `--verbose` / `-v` CLI flag

In `core/cli.py`, add:

```python
parser.add_argument(
    "-v", "--verbose",
    action="store_true",
    default=False,
    help="Enable verbose output with debug-level details",
)
```

In `main.py`, call `setup_logging(verbose=args.verbose)` before any other logic.

### 3.5 Final result output

Only the result line (`CLANS file: ...`) should remain as plain `print()` — it is the program's primary output and should always appear, regardless of log configuration. Format them distinctly:

```python
print(f"  CLANS file:    {clans_file_path}")
```

### 3.6 Progress indicators (Issue #30)

For long-running operations, add progress counters at `INFO` level:

**Structure downloads** (`structure_utils.py`):  
Print a running counter:

```
[INFO   ] Downloading structures... (50 proteins)
[INFO   ] Downloaded 10/50 structures...
[INFO   ] Downloaded 20/50 structures...
...
[INFO   ] Downloaded 48/50 structures (2 failed)
```

Implementation options:

- **Option A (simple):** Log every N-th download (e.g., every 10 or every 20%)
- **Option B (rich):** Use a progress bar library like `tqdm` (new dependency) behind a wrapper

**Similarity computation** (`struct_sim_computer.py`):  
Since tool execution is a single subprocess call, we can only report start/end. The tool's own stderr could be passed through at `DEBUG` level.

**CLANS file generation** (`clans_file_generator.py`):  
Fast operation — just start/end messages are sufficient.

---

## 4. Implementation Plan

### Step 1: Create `src/clans3d/utils/log.py`

- Define `setup_logging(verbose: bool)` function
- Consistent `[LEVEL] message` format

### Step 2: Add `-v` / `--verbose` flag to CLI

- Update `core/cli.py` with the new argument
- Call `setup_logging()` at the top of `main()`

### Step 3: Replace `print()` with `logger` calls module by module

Work through each module, replacing prints with the appropriate log level. Order:

1. **`main.py`** — setup logging, keep final result as `print()`
2. **`core/clans_file_generator.py`** — `logger.info()` for step start/complete
3. **`utils/dependency_checks.py`** — `logger.debug()` for tool path info
4. **`utils/file_utils.py`** — `logger.warning()` for download failures
5. **`utils/structure_utils.py`** — `logger.info()` for progress, `logger.warning()` for failures
6. **`utils/fasta_utils.py`** — `logger.warning()` for missing sequences, `logger.debug()` for filtering steps
7. **`similarity/struct_sim_computer.py`** — `logger.info()` for step flow, `logger.warning()` for score mismatches
8. **`similarity/struct_sim_tool.py`** — `logger.error()` for subprocess failures
9. **`similarity/foldseek.py`** — `logger.error()` for tool errors, `logger.debug()` for stdout/stderr
10. **`similarity/tmalign.py`** — `logger.error()` for parse failures
11. **`similarity/usalign.py`** — `logger.error()` for parse failures
12. **`legacy/utils_old_clans.py`** — `logger.debug()` for command args
13. **`benchmark/benchmark.py`** — Migrate to logger but keep structured table formatting

### Step 4: Add progress reporting for structure downloads

- Add download counter in `process_fasta_file()` and `process_tsv_file()`
- Log at intervals (every 10 items or every 20%, whichever is smaller)

### Step 5: Propagate `verbose` to library usage (non-CLI)

- `PipelineConfig` gets an optional `verbose: bool = False` field
- `ClansPipeline.__init__` calls `setup_logging()` if not already configured
- `Benchmark` can set its own verbosity level

### Step 6: Testing & review

- Verify default output is clean and ordered for a typical run
- Verify `--verbose` adds debug details without breaking format
- Verify benchmark output still renders its tables correctly
- Ensure no bare `print()` remains in the codebase (except final result in `main.py` and benchmark tables)

---

## 5. Files Changed Summary

| File                                            | Action                                             |
| ----------------------------------------------- | -------------------------------------------------- |
| `src/clans3d/utils/log.py`                      | **New** — logging setup                            |
| `src/clans3d/core/cli.py`                       | Add `-v` flag                                      |
| `src/clans3d/main.py`                           | Init logging, keep result prints                   |
| `src/clans3d/core/pipeline.py`                  | Accept `verbose` in config                         |
| `src/clans3d/core/clans_file_generator.py`      | `print` → `logger.info`                            |
| `src/clans3d/utils/dependency_checks.py`        | `print` → `logger.debug`                           |
| `src/clans3d/utils/file_utils.py`               | `print` → `logger.warning`                         |
| `src/clans3d/utils/structure_utils.py`          | `print` → `logger.info/warning` + progress counter |
| `src/clans3d/utils/fasta_utils.py`              | `print` → `logger.warning/debug`                   |
| `src/clans3d/similarity/struct_sim_computer.py` | `print` → `logger.info/warning`                    |
| `src/clans3d/similarity/struct_sim_tool.py`     | `print` → `logger.error`                           |
| `src/clans3d/similarity/foldseek.py`            | `print` → `logger.error/debug`                     |
| `src/clans3d/similarity/tmalign.py`             | `print` → `logger.error`                           |
| `src/clans3d/similarity/usalign.py`             | `print` → `logger.error`                           |
| `src/clans3d/legacy/utils_old_clans.py`         | `print` → `logger.debug`                           |

---

## 6. Decisions

1. **Progress bars:** Log-line counters (no `tqdm` dependency).
2. **Benchmark module:** Keeps its own `print()` tables — no migration to logger.
3. **Final output:** Only the CLANS file path is printed as a plain `print()` line.
4. **Quiet mode:** Deferred to a follow-up; the verbose/default split covers the immediate need.
