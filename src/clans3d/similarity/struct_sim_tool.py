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
    Base class for structural similarity tools.
    
    Subclasses must implement :meth:`start_run`, :meth:`_parse_output`,
    and :meth:`_log_progress`.
    """
    def __init__(self, name: str, description: str, working_dir: str):
        self.name = name
        self.description = description
        self.working_dir = working_dir
        self.command = []
        self.output = None
        self.expected_number_of_scores = 0
        
        
    def start_run(self, structures_dir: str, expected_number_of_scores: int) -> pd.DataFrame:
        """
        Initializes the self.command list with the necessary parameters to run the tool and then returns _execute_run.
        This method should be overridden by subclasses to implement specific tool logic.
        """
        raise NotImplementedError("Subclasses should implement this method to run the tool.")
    
    
    def _execute_run(self, expected_number_of_scores: int) -> pd.DataFrame:
        """
        Executes ``self.command`` via ``Popen`` and streams the output,
        periodically calling :meth:`_log_progress` so the user can
        track that the tool is still running.
        
        Args:
            expected_number_of_scores: The expected number of pairwise
                scores.  Stored on the instance and available to
                :meth:`_log_progress` for percentage calculation.
        
        Returns:
            A DataFrame with the parsed similarity scores, or an empty
            DataFrame on error.
        """
        self.expected_number_of_scores = expected_number_of_scores
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
        Block until *process* finishes, calling :meth:`_log_progress`
        every ``PROGRESS_LOG_INTERVAL`` seconds.

        A final call is made after the process exits so the user always
        sees the last status.

        Args:
            process:        The running subprocess to monitor.
            stdout_reader:  Reader dict for the stdout pipe.
            stderr_reader:  Reader dict for the stderr pipe.
        """
        while process.poll() is None:
            time.sleep(PROGRESS_LOG_INTERVAL)
            self._log_progress(stdout_reader, stderr_reader)
        # One final log after the process exits
        self._log_progress(stdout_reader, stderr_reader)


    def _log_progress(self, stdout_reader: dict,
                      stderr_reader: dict) -> None:
        """
        Log a progress message based on the current subprocess output.

        Args:
            stdout_reader: Reader dict for the stdout pipe.
            stderr_reader: Reader dict for the stderr pipe.
        """
        raise NotImplementedError("Subclasses should implement this method to log progress during execution.")

    
    def _parse_output(self) -> pd.DataFrame:
        """
        Parses the output of the tool to extract the similarity scores.
        This method should be overridden by subclasses to implement specific parsing logic.
        """
        raise NotImplementedError("Subclasses should implement this method to parse output.")
