"""
CLI argument parsing for Clans-3D.

Provides the argparse setup for the main entry point and the config-file
pre-parser.
"""

import sys
import argparse

from clans3d.core.input_file_type import InputFileType
from clans3d.core.config_file import ConfigFile
from clans3d.similarity.tool_type import ToolType


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
        choices=[t.value for t in InputFileType],
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
        "-v", "--verbose",
        action="store_true",
        default=False,
        help="enable verbose output with debug-level details",
    )
    
    parser.add_argument(
        "-q", "--quiet",
        action="store_true",
        default=False,
        help="disable all output except for errors",
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
