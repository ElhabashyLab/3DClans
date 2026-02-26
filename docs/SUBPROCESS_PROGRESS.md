# Subprocess Progress Reporting — Implementation Plan

**Goal:** Show the user progress (in 10% increments) while external similarity tools (Foldseek, USalign) run as subprocesses.

---

## Problem

`StructSimTool._execute_run()` uses `subprocess.run(capture_output=True)`, which **buffers all output** until the process finishes. The user sees nothing until the entire computation is done — problematic for large datasets where runs take minutes to hours.

---

## Approach

### 1. Replace `subprocess.run()` with `subprocess.Popen()` in `_execute_run()`

Switch to `Popen` so we can **stream stderr in real-time** while still capturing stdout.

**File:** `src/clans3d/similarity/struct_sim_tool.py`

```python
# Current (blocking, no progress):
result = subprocess.run(self.command, capture_output=True, text=True, check=True)
self.output = result.stdout

# New (streaming):
process = subprocess.Popen(
    self.command,
    stdout=subprocess.PIPE,
    stderr=subprocess.PIPE,
    text=True,
)
# Stream stderr for progress while accumulating stdout
self.output, _ = self._stream_with_progress(process)
```

### 2. Add `_stream_with_progress()` to `StructSimTool`

Reads stderr line-by-line (split on `\r` and `\n`), extracts percentage values via regex, and logs at every 10% milestone.

```python
def _stream_with_progress(self, process):
    """
    Reads stdout/stderr from a Popen process.
    Parses percentage from stderr and logs progress at 10% intervals.
    Returns (stdout_content, stderr_content).
    """
```

**Key details:**

- Use a **background thread** to read stdout into a buffer (prevents deadlock when both pipes are used)
- Read stderr character-by-character in the main thread, splitting on `\r`/`\n` to detect complete progress lines
- Parse percentage via regex: `r'(\d+\.\d+)%'`
- Track last reported milestone (0, 10, 20, …, 100) and only log when a new 10% threshold is crossed
- Check process return code after completion, raise on non-zero

### 3. Add `_on_progress()` hook in `StructSimTool`

```python
def _on_progress(self, phase: str, percent: float):
    """Called when a new 10% progress milestone is reached."""
    logger.info("%s – %s: %d%%", self.name, phase, int(percent))
```

Subclasses can override this if they need custom behavior.

### 4. Detect phases from Foldseek stderr

Foldseek emits labeled text lines between progress bars. Parse these to set the current phase:

| Stderr text (non-progress line)            | Phase label         |
| ------------------------------------------ | ------------------- |
| `Index table: counting k-mers`             | indexing            |
| `Starting prefiltering scores calculation` | prefiltering        |
| `structurealign`                           | structure alignment |

**Implementation:** In `_stream_with_progress()`, when a stderr line doesn't contain `[` or `%`, check it against known phase keywords and update the current phase label.

Example output the user will see:

```
[INFO   ] foldseek – prefiltering: 10%
[INFO   ] foldseek – prefiltering: 20%
...
[INFO   ] foldseek – structure alignment: 50%
[INFO   ] foldseek – structure alignment: 60%
...
[INFO   ] foldseek – structure alignment: 100%
```

### 5. USalign progress via stdout line counting

USalign writes one result line per pair to stdout. We know the expected pair count from the number of structure files (`n*(n-1)/2`).

**File:** `src/clans3d/similarity/usalign.py`

- Override `_execute_run()` (or use a flag) so that stdout is read line-by-line
- Pass `expected_pairs` to the method
- Count lines, compute percentage, log at 10% milestones
- Accumulate all stdout lines into `self.output` for subsequent parsing

```python
def start_run(self, structures_dir):
    ...
    n = len(structure_files)
    self._expected_pairs = n * (n - 1) // 2
    return self._execute_run()
```

Example output:

```
[INFO   ] USalign – aligning pairs: 10%
[INFO   ] USalign – aligning pairs: 50%
[INFO   ] USalign – aligning pairs: 100%
```

### 6. Foldseek helper subprocesses (`_create_database`, `_parse_output`)

**File:** `src/clans3d/similarity/foldseek.py`

These two methods also call `subprocess.run()`:

- `_create_database()` — usually fast, but log a single start/done message
- `_parse_output()` (convertalis) — also usually fast, same treatment

For consistency, either:

- **(a)** Reuse `_stream_with_progress()` for these too, or
- **(b)** Keep `subprocess.run()` but add `logger.info()` messages before/after (simpler, since these steps are typically < 1 second)

**Recommendation:** Option **(b)** — keep it simple for fast steps.

```python
logger.info("Creating Foldseek database...")
# subprocess.run(...)
logger.info("Foldseek database created.")
```

---

## Files to Change

| File                                        | Changes                                                                                    |
| ------------------------------------------- | ------------------------------------------------------------------------------------------ |
| `src/clans3d/similarity/struct_sim_tool.py` | Add `_stream_with_progress()`, `_on_progress()`; rewrite `_execute_run()` to use `Popen`   |
| `src/clans3d/similarity/foldseek.py`        | Add phase detection keywords; add info logs around `_create_database` and convertalis      |
| `src/clans3d/similarity/usalign.py`         | Override `_execute_run()` or adapt `start_run()` to count stdout lines and report progress |

No changes needed in `struct_sim_computer.py`, `main.py`, or `log.py`.

---

## No New Dependencies

Everything is implemented with:

- `subprocess.Popen` (stdlib)
- `threading.Thread` (stdlib, for reading stdout in background)
- `re` (stdlib, for parsing `%`)
- Existing `logging` module

Quiet mode (`-q`) automatically suppresses progress since `logger.info()` is silenced at `ERROR` level.

---

## Testing

1. **Manual test with 33 structures** (already in `work/structures/`):

   ```bash
   clans3d -l examples/small_fasta_files/50.fasta -i fasta -t foldseek
   ```

   Verify 10% milestone logs appear during the search phase.

2. **Manual test with USalign:**

   ```bash
   clans3d -l examples/small_fasta_files/5.fasta -i fasta -t usalign
   ```

3. **Edge cases:**
   - Very fast runs (< 1 second): may jump from 0% to 100% — acceptable
   - Tool not found: existing error handling unchanged
   - Quiet mode: no progress output — verify with `-q`
