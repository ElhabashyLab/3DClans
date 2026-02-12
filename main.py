import sys
import argparse
import os
from structure_utils import fetch_pdbs
from fasta_utils import generate_fasta_from_uids_with_regions
from StructSimComputer import StructSimComputer
from ToolType import ToolType
from ClansFileGenerator import ClansFileGenerator
from InputFileType import InputFileType
from ConfigFile import ConfigFile
from dependency_checks import verify_tool_dependencies


def _set_up_parser() -> argparse.ArgumentParser:
    """
    Sets up the command line parser.
    :return: the parser
    """
    parser = argparse.ArgumentParser(description="Command line parser to read fasta A2M or tsv files and generating clans files based on structural similarity scores.")
    
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
        "-c", "--conf",
        required=False,
        type=str,
        help="specifies the path to the configuration file with the format -<key> <value>. The configuration file arguments will be overwritten by command line arguments."
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
        default=None,
        help="specifies the scoring method to use for Foldseek (default: evalue)",
    )
    return parser


def _set_up_conf_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Parser for configuration file.")
    parser.add_argument("-c", "--conf", type=str)
    return parser

    
def main():
    # parse configuration file arguments and merge with command line arguments (command line arguments take precedence)
    conf_parser = _set_up_conf_parser()
    conf_args, _ = conf_parser.parse_known_args()
    config_argv = []
    if conf_args.conf:
        conf_dict = ConfigFile(conf_args.conf).read_config()
        config_argv = _config_dict_to_argv(conf_dict)    
    merged_argv = config_argv + sys.argv[1:]
    # parse merged arguments
    parser = _set_up_parser()
    args = parser.parse_args(merged_argv)
    # check dependencies for the selected tool
    verify_tool_dependencies(ToolType(args.tool))
    # process input file and create clans file
    input_file = args.load
    clans_file_path, cleaned_input_file_path = create_clans_file(
        input_file,
        InputFileType(args.input_type),
        ToolType(args.tool),
        args.score,
    )
    return clans_file_path, cleaned_input_file_path
    
    
def _config_dict_to_argv(config: dict) -> list[str]:
    """
    Converts a config dictionary to a list of argv strings.
    The arguments are formatted as '--key value'.

    Args:
        config (dict): A dictionary containing configuration key-value pairs.

    Returns:
        list[str]: A list of strings formatted as command-line arguments.
    """
    argv = []
    for key, value in config.items():
        argv.append(f"--{key}")
        argv.append(str(value))
    return argv
    
    
def create_clans_file(
    input_file_path: str,
    input_file_type: InputFileType,
    selected_tool: ToolType,
    foldseek_score: str | None,
    structures_dir: str = "structures",
    out_dir_path: str = "clans_files"
    ) -> tuple[str, str]:
    """
    Creates a clans file.

    Args:
        input_file_path (str): path to the saved input file.
        input_file_type (InputFileType): Type of the input file.
        selected_tool (ToolType): Type of tool which is to generate structural similarity scores.
        foldseek_score (str): Specifies which score to use if foldseek is the selected tool. Otherwise is None.
        structures_dir (str): Path to the directory containing the structures.
        out_dir_path (str): Path to the output directory where the clans file will be saved.

    Returns:
        (str, str): Path to generated clans file and path to cleaned input file as fasta.
    """
    scores_computer = _set_up_scores_computer(selected_tool, foldseek_score)
    uids_with_regions = fetch_pdbs(input_file_path, input_file_type, structures_dir)
    scores = scores_computer.run(selected_tool, structures_dir)
    
    input_file_name = os.path.basename(input_file_path).split(".")[0]
    os.makedirs("input_file_storage", exist_ok=True)
    cleaned_input_file_path = os.path.join("input_file_storage", f"{input_file_name}_cleaned.fasta")
    if input_file_type == InputFileType.FASTA:
        cleaned_input_file_path = generate_fasta_from_uids_with_regions(uids_with_regions, cleaned_input_file_path, input_file_path)
    else:
        cleaned_input_file_path = generate_fasta_from_uids_with_regions(uids_with_regions, cleaned_input_file_path)
    
    clans_generator = ClansFileGenerator()
    input_file_cleaned_name = os.path.basename(cleaned_input_file_path).split(".")[0]
    out_path = os.path.join(out_dir_path, f"{input_file_cleaned_name}.clans")
    clans_file_path = clans_generator.generate_clans_file(scores, cleaned_input_file_path, out_path)
    return clans_file_path, cleaned_input_file_path
    

def _set_up_scores_computer(selected_tool: ToolType, foldseek_score: str | None) -> StructSimComputer:
    """Creates an instance of a StructSimComputer with the given selected_tool, and a possible foldseekscore.
    
    Args:
        selected_tool (ToolType): The ToolType used by the StructSimComputer
        foldseek_score (str): The foldseek_score used by the StrucSimComputer. Is None if selected_tool is not foldseek

    Raises:
        ValueError: Is raised if a foldseek_score is specified but foldseek not selected_tool

    Returns:
        StructSimComputer: The instance of the StructSimComputer
    """
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
