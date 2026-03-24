"""
CLI argument parsing for 3DClans.

Provides the argparse setup for the main entry point and the config-file
pre-parser.
"""

import os
import sys
import argparse

from clans3d.core.input_file_type import InputFileType
from clans3d.core.config_file import ConfigFile
from clans3d.similarity.tool_type import ToolType
from clans3d.similarity.tm_mode import TmMode

_DEFAULT_OUTPUT_DIR = os.path.join("output", "clans_files")


def _build_main_parser() -> argparse.ArgumentParser:
    """Build the full CLI argument parser."""
    parser = argparse.ArgumentParser(
        description=(
            "Command line parser to read fasta A2M or tsv files and "
            "generating clans files based on structural similarity scores."
        )
    )

    parser.add_argument(
        "-l", "--load",
        required=True,
        type=str,
        help="specifies the path to the input file",
    )

    parser.add_argument(
        "-i", "--input_type",
        required=True,
        choices=[InputFileType.FASTA.value, InputFileType.A2M.value, 
                 InputFileType.A3M.value, InputFileType.TSV.value],
        help="specifies the type of the input file",
    )

    parser.add_argument(
        "-c", "--conf",
        required=False,
        type=str,
        help=(
            "specifies the path to the configuration file with the format "
            "-<key> <value>. The configuration file arguments will be "
            "overwritten by command line arguments."
        ),
    )

    parser.add_argument(
        "-t", "--tool",
        required=True,
        choices=[t.value for t in ToolType],
        help="specifies the tool to use for computing similarity scores",
    )

    parser.add_argument(
        "-s", "--score",
        required=False,
        choices=["evalue", "TM"],
        default=None,
        help="specifies the scoring method to use for Foldseek (default: evalue)",
    )

    parser.add_argument(
        "-m", "--tm_mode",
        required=False,
        choices=[m.value for m in TmMode],
        default=TmMode.MIN.value,
        help=(
            "specifies how qtm and ttm are combined into the final distance score: "
            "min, max, or mean (default: min). "
            "Used by USalign and Foldseek when -s TM."
        ),
    )

    verbosity = parser.add_mutually_exclusive_group()
    verbosity.add_argument(
        "-v", "--verbose",
        required=False,
        action="store_true",
        default=False,
        help="enable verbose output with debug-level details",
    )

    verbosity.add_argument(
        "-q", "--quiet",
        required=False,
        action="store_true",
        default=False,
        help="disable all output except for errors",
    )

    parser.add_argument(
        "-o", "--out",
        required=False,
        type=str,
        default=None,
        help=(
            "output path for the generated CLANS file. "
            "Can be a full file path (e.g. results/my_output.clans) "
            "or a directory (e.g. results/). "
            "When omitted the file is written to output/clans_files/."
        ),
    )
    
    def _positive_int(value):
        ivalue = int(value)
        if ivalue <= 0:
            raise argparse.ArgumentTypeError(f"{value} is not a positive integer")
        return ivalue

    parser.add_argument(
        "-w", "--workers",
        required=False,
        type=_positive_int,
        default=10,
        metavar="N",
        help=(
            "number of parallel threads for structure downloads (default: 10). "
            "Increase for faster downloads on large datasets; "
            "decrease if the AlphaFold API starts rate-limiting."
        ),
    )

    return parser

def _build_conf_parser() -> argparse.ArgumentParser:
    """Build a pre-parser to detect the config file argument."""
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument(
        "-c", "--conf",
        type=str,
        help="specifies the path to the configuration file with the format -<key> <value>",
    )
    return parser

def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    """Parse CLI arguments, merging an optional config file.

    If ``-c / --conf`` is present the config file is read first and its
    entries are placed *before* the real CLI arguments so that the latter
    take precedence.

    Args:
        argv: Argument list to parse.  Defaults to ``sys.argv[1:]``.

    Returns:
        Parsed :class:`argparse.Namespace`.
    """
    # grab argv from cli
    if argv is None:
        argv = sys.argv[1:]

    # Pre-parse to detect a config file
    conf_parser = _build_conf_parser()
    conf_args, _ = conf_parser.parse_known_args(argv)
    path_to_config = conf_args.conf

    config_argv: list[str] = []
    if conf_args.conf:
        config_argv = ConfigFile(path_to_config).config_to_argv()

    # Merge config arguments with CLI arguments, giving precedence to CLI
    merged_argv = config_argv + argv

    parser = _build_main_parser()
    return parser.parse_args(merged_argv)


def resolve_output_path(out: str | None) -> tuple[str, str | None]:
    """Resolve the ``-o / --out`` argument into output_dir and output_filename.

    Args:
        out: The raw value of ``args.out`` (may be ``None``).

    Returns:
        Tuple of ``(output_dir, output_filename)``.  *output_filename* is
        ``None`` when only a directory was given or *out* was ``None``.
    """
    if out is None:
        return _DEFAULT_OUTPUT_DIR, None

    basename = os.path.basename(out)
    _, ext = os.path.splitext(basename)
    if ext == ".clans":
        return os.path.dirname(out) or ".", basename
    return out, None
