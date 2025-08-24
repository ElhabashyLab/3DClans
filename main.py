import sys
import argparse
import os
from fasta2PDB import *
from StructSimComputer import StructSimComputer
from ToolType import ToolType
from ClansFileGenerator import ClansFileGenerator


def _save_file(file_path, is_clans):
    """
    Reads and stores a CLANS or FASTA file with a given file-path.
    :param file_path: path to file
    :param is_clans: specifies if input is CLANS or FASTA file
    """
    print(f"processing of {file_path}:")
    try:
        with open(file_path, 'r') as file:
            content = file.read()
        os.makedirs('input_file_storage', exist_ok=True)
        if is_clans:
            with open('input_file_storage/input_file.clans', 'w') as input_file:
                input_file.write(content)
        else:
            with open('input_file_storage/input_file.fasta', 'w') as input_file:
                input_file.write(content)
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.", file=sys.stderr)
    except IOError as e:
        print(f"Error reading file '{file_path}': {e}", file=sys.stderr)


def _set_up_parser():
    """
    Sets up the command line parser.
    :return: the parser
    """
    parser = argparse.ArgumentParser(description="command line parser to read fasta or clans file")
    file_type = parser.add_mutually_exclusive_group(required=True)
    
    file_type.add_argument(
        "--f", "-fasta",
        type=str,
        help="specifies the input as fasta file"
    )
    
    file_type.add_argument(
        "--c", "-clans",
        type=str,
        help="specifies the input as clans file"
    )
    
    parser.add_argument(
        "--t", "-tool",
        required=True,
        choices=[tool.value for tool in ToolType],
        help="specifies the tool to use for computing similarity scores"
    )
    return parser


def main():
    parser = _set_up_parser()
    # reading arguments
    args = parser.parse_args()
    selected_tool = ToolType(args.t)
    if args.f:
        _save_file(args.f, False)
        # fasta to pdb conversion
        input_file = "input_file_storage/input_file.fasta"
        cleaned_input_file = fetch_pdbs_from_uids(input_file, "PDBs")
        # pdb to scores conversion
        computer = StructSimComputer()
        scores = computer.run(selected_tool, "PDBs")
        # clans file generation
        generator = ClansFileGenerator()
        clans_file_path = generator.generate_clans_file(scores, cleaned_input_file)
    else:
        _save_file(args.c, True)
        raise NotImplementedError

    
if __name__ == "__main__":
    main()
