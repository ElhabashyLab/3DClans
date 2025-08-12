import sys
import argparse
import os

def read_file(file_path, is_clans):
    """
    Reads and stores a CLANS or FASTA file with a given file-path.
    :param file_path: path to file
    :param is_clans: specifies if input is CLANS or FASTA file
    """
    print(f"processing of {file_path}:")
    try:
        with open(file_path, 'r') as file:
            content = file.read()
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


def main():
    parser = argparse.ArgumentParser(description="command line parser to read fasta or clans file")
    file_type = parser.add_mutually_exclusive_group(required=True)

    # arguments for file_type
    file_type.add_argument(
        "--f", "-fasta",
        action="store_true",
        help="specifies the input as fasta file"
    )
    file_type.add_argument(
        "--c", "-clans",
        action="store_true",
        help="specifies the input as clans file"
    )

    # arguments for parser
    parser.add_argument(
        "file",
        type=str,
        help="represents path to input file"
    )

    # reading arguments
    args = parser.parse_args()
    if args.f:
        read_file(args.file, False)
    else:
        read_file(args.file, True)

    # fasta to pdb conversion
    # pdb to scores conversion
    # clans file generation
    # saving clans file
if __name__ == "__main__":
    main()
