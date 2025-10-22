import sys
import argparse
import os
from utils_for_PDB import *
from StructSimComputer import StructSimComputer
from ToolType import ToolType
from ClansFileGenerator import ClansFileGenerator
from InputFileType import InputFileType
    

def _save_file(file_path: str, input_file_type: InputFileType) -> str:
    """
    Reads and stores a CLANS file, FASTA file or a list of uids with a given file-path.
    
    :param file_path: path to file
    :param input_file_type: specifies if input is CLANS file, FASTA file or a list of uids as a tsv file
    :return: path to saved file
    :raises: FileNotFoundError, IOError if file operations fail
    """
    print(f"Processing {file_path}...")
    try:
        with open(file_path, 'r') as file:
            content = file.read()
    except FileNotFoundError:
        print(f"Error: File '{file_path}' not found.", file=sys.stderr)
        raise FileNotFoundError
    except IOError as e:
        print(f"Error reading file '{file_path}': {e}", file=sys.stderr)
        raise IOError
    
    os.makedirs('input_file_storage', exist_ok=True)
    output_filename = f"input_file.{input_file_type.value}"
    path_to_stored_file = os.path.join('input_file_storage', output_filename)
    
    try:
        with open(path_to_stored_file, 'w') as input_file:
            input_file.write(content)
        print(f"File successfully saved to: {path_to_stored_file}")
        return path_to_stored_file
    except IOError as e:
        print(f"Error writing file '{path_to_stored_file}': {e}", file=sys.stderr)
        raise IOError


def _set_up_parser() -> argparse.ArgumentParser:
    """
    Sets up the command line parser.
    :return: the parser
    """
    parser = argparse.ArgumentParser(description="command line parser to read fasta or clans file")
    
    parser.add_argument(
        "-l", "--load",
        required=True,
        type=str,
        help="specifies the path to the input file"
    )
    
    parser.add_argument(
        "-i", "--input_type",
        required=True,
        choices=[type.value for type in InputFileType],
        help="specifies the type of the input file"
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
    """
    The main method.
    Reads the input of the user and directs the generation of a clans file
    """
    parser = _set_up_parser()
    args = parser.parse_args()
    input_file_type = InputFileType(args.input_type)
    path_to_input_file = args.load
    selected_tool = ToolType(args.tool)
    foldseek_score = args.score
    saved_input_file = _save_file(path_to_input_file, input_file_type)
    clans_file = _create_clans_file(saved_input_file, input_file_type, selected_tool, foldseek_score)
            

def _create_clans_file(saved_input_file, input_file_type, selected_tool, foldseek_score):
    """
    Creates a clans file.

    Args:
        saved_input_file (str): path to the saved input file.
        input_file_type (InputFileType): Type of the input file.
        selected_tool (ToolType): Type of tool which is to generate structural similarity scores.
        foldseek_score (str): Specifies which score to use if foledseek is the selected tool. Otherwise is None.

    Returns:
        str: Path to generated clans file.
    """
    cleaned_input_file = fetch_pdbs(saved_input_file, input_file_type, "PDBs")
    scores_computer = _set_up_scores_computer(selected_tool, foldseek_score)
    scores = scores_computer.run(selected_tool, "PDBs")
    clans_generator = ClansFileGenerator()
    clans_file_path = clans_generator.generate_clans_file(scores, cleaned_input_file)
    return clans_file_path
    

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
