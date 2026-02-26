import logging
import subprocess
import threading
import time
import pandas as pd

logger = logging.getLogger(__name__)

# Interval in seconds between progress log messages during subprocess execution
PROGRESS_LOG_INTERVAL = 5


class StructSimTool():
    """
    This module defines a base class for tools that can be used in a benchmarking context.
    """
    def __init__(self, name: str, description: str, working_dir: str):
        self.name = name
        self.description = description
        self.working_dir = working_dir
        self.command = []
        self.output = None
        
        
    def start_run(self, structures_dir: str) -> pd.DataFrame:
        """
        Initializes the self.command list with the necessary parameters to run the tool and then returns _execute_run.
        This method should be overridden by subclasses to implement specific tool logic.
        """
        raise NotImplementedError("Subclasses should implement this method to run the tool.")
    
    
    def _execute_run(self) -> pd.DataFrame:
        """
        Executes self.command in a subprocess and captures the output.
        Streams stderr and logs it at regular intervals so the user
        can see that the tool is still running.
        
        Returns:
            pd.DataFrame: A DataFrame containing the similarity scores between the structures.
        """
        logger.debug("Running command: %s", " ".join(self.command))
        try:
            process = subprocess.Popen(
                self.command,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )
            self.output, stderr = self._stream_with_progress(process)
            if process.returncode != 0:
                raise subprocess.CalledProcessError(
                    process.returncode, self.command,
                    output=self.output, stderr=stderr,
                )
            scores = self._parse_output()
        except subprocess.CalledProcessError as e:
            logger.error("Error running %s with %s: %s", self.name, self.command, e)
            scores = pd.DataFrame()  # Return an empty DataFrame on error
        return scores


    def _stream_with_progress(self, process: subprocess.Popen) -> tuple[str, str]:
        """
        Reads stdout and stderr from a running ``Popen`` process and
        periodically logs the latest output line so the user can track
        that the subprocess is still active.

        Both pipes are consumed by background reader threads to avoid
        deadlocks.  The main thread polls the process at
        ``PROGRESS_LOG_INTERVAL``-second intervals and logs the most
        recent non-empty line from either stream.

        Args:
            process: The running subprocess whose output should be
                streamed and logged.

        Returns:
            A ``(stdout, stderr)`` tuple with the full captured output.
        """
        stdout_reader = self._start_pipe_reader(process.stdout)
        stderr_reader = self._start_pipe_reader(process.stderr)

        self._poll_and_log(process, stdout_reader, stderr_reader)

        stdout_reader["thread"].join()
        stderr_reader["thread"].join()
        process.wait()

        return (
            "".join(stdout_reader["chunks"]),
            "".join(stderr_reader["chunks"]),
        )


    def _start_pipe_reader(self, pipe) -> dict:
        """
        Start a daemon thread that reads lines from *pipe* into a
        shared dictionary.

        The returned dict contains:

        - ``chunks``  – list of raw lines read so far.
        - ``latest_line`` – the most recent non-blank line (stripped).
        - ``thread`` – the ``Thread`` object (already started).

        Args:
            pipe: A file-like pipe object (e.g. ``process.stdout``).

        Returns:
            A dict with the accumulated output and the reader thread.
        """
        reader: dict = {"chunks": [], "latest_line": "", "thread": None}

        def _read():
            for line in pipe:
                reader["chunks"].append(line)
                stripped = line.rstrip()
                if stripped:
                    reader["latest_line"] = stripped
            pipe.close()

        thread = threading.Thread(target=_read, daemon=True)
        reader["thread"] = thread
        thread.start()
        return reader


    def _poll_and_log(self, process: subprocess.Popen,
                      stdout_reader: dict, stderr_reader: dict) -> None:
        """
        Block until *process* finishes, logging progress at regular
        intervals.

        Every ``PROGRESS_LOG_INTERVAL`` seconds the latest captured
        output line is emitted via :func:`logging.info`.  ``stderr`` is
        preferred (where tools like Foldseek write progress bars), with
        a fallback to ``stdout`` (used by USalign).

        A final log entry is emitted once the process exits so the user
        always sees the last status.

        Args:
            process:        The running subprocess to monitor.
            stdout_reader:  Reader dict returned by
                :meth:`_start_pipe_reader` for stdout.
            stderr_reader:  Reader dict returned by
                :meth:`_start_pipe_reader` for stderr.
        """
        while process.poll() is None:
            time.sleep(PROGRESS_LOG_INTERVAL)
            self._log_latest(stdout_reader, stderr_reader)
        # One final log after the process exits
        self._log_latest(stdout_reader, stderr_reader)


    def _log_latest(self, stdout_reader: dict,
                    stderr_reader: dict) -> None:
        """
        Log the most recent non-empty output line.

        Prefers the ``stderr`` stream (common for progress output),
        falling back to ``stdout`` if stderr has nothing new.

        Args:
            stdout_reader: Reader dict for the stdout pipe.
            stderr_reader: Reader dict for the stderr pipe.
        """
        latest = stderr_reader["latest_line"] or stdout_reader["latest_line"]
        if latest:
            logger.info("%s: %s", self.name, latest)

    
    def _parse_output(self) -> pd.DataFrame:
        """
        Parses the output of the tool to extract the similarity scores.
        This method should be overridden by subclasses to implement specific parsing logic.
        """
        raise NotImplementedError("Subclasses should implement this method to parse output.")
