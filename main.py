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
    :return: path to saved file
    """
    print(f"processing of {file_path}:")
    try:
        with open(file_path, 'r') as file:
            content = file.read()
        os.makedirs('input_file_storage', exist_ok=True)
        if is_clans:
            path_to_stored_file = 'input_file_storage/input_file.clans'
            with open(path_to_stored_file, 'w') as input_file:
                input_file.write(content)
        else:
            path_to_stored_file = 'input_file_storage/input_file.fasta'
            with open(path_to_stored_file, 'w') as input_file:
                input_file.write(content)
        return path_to_stored_file
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
        "-f", "--fasta",
        type=str,
        help="specifies the input as fasta file"
    )
    
    file_type.add_argument(
        "-c", "--clans",
        type=str,
        help="specifies the input as clans file"
    )
    
    parser.add_argument(
        "-t", "--tool",
        required=True,
        choices=[tool.value for tool in ToolType],
        help="specifies the tool to use for computing similarity scores"
    )
    
    parser.add_argument(
        "-s", "--score",
        required=False,
        choices=["evalue", "TM"],
        default="evalue",
        help="specifies the scoring method to use for Foldseek (default: evalue)",
        
    )
    return parser


def main():
    parser = _set_up_parser()
    args = parser.parse_args()
    # process input as clans file
    if args.clans:
        input_file = _save_file(args.clans, True)
        raise NotImplementedError
    # process input as fasta file
    elif args.fasta:
        input_file = _save_file(args.fasta, False)
        selected_tool = ToolType(args.tool)
        foldseek_score = args.score
        scores_computer = _set_up_scores_computer(selected_tool, foldseek_score)
        clans_generator = ClansFileGenerator()
        cleaned_input_file = fetch_pdbs(input_file, "PDBs")
        scores = scores_computer.run(selected_tool, "PDBs")
        clans_file_path = clans_generator.generate_clans_file(scores, cleaned_input_file)
    else:
        raise ValueError("Please specify either a fasta or clans file as input.")
            

def _set_up_scores_computer(selected_tool, foldseek_score):
    if selected_tool == ToolType.FOLDSEEK:
        if foldseek_score == None or foldseek_score == "evalue":
            return StructSimComputer()
        else:
            return StructSimComputer("TM")
    elif selected_tool != ToolType.FOLDSEEK and foldseek_score != None:
        raise ValueError("The foldseek score (-s --score) can only be specified if foldseek is chosen as tool.")
    else:
        return StructSimComputer()


if __name__ == "__main__":
    main()
