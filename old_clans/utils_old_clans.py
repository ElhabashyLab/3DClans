import subprocess


# This is the path to the recovered clans project directory (github-repo: https://github.com/AronWichtner/clans-recovered.git)


def run_clans_headless(recovered_clans_path: str, input_file: str, output_file: str, rounds: int = 100):
    """
    Runs recovered clans in headless mode on the given clans file and saves the output to the output file.
    Args:
        input_file: A path to the input clans file
        output_file: A path to the output clans file
        recovered_clans_path: A path to the recovered clans project directory
    Returns: None
    """
    gradlew = "./gradlew"
    args = [
        gradlew, "run", "--no-daemon",
        f'--args=-nographics T -load {input_file} -dorounds {rounds} -saveto {output_file}'
    ]
    subprocess.run(args, cwd=recovered_clans_path, check=True)


def run_multiple_clans_headless(recovered_clans_path: str, input_output_files: dict, rounds: int = 100):
    """
    Runs recovered clans in headless mode on each of the given clans files and saves the output to the corresponding output files.
    The input_output_files dict should contain input_file: output_file pairs.
    Args:
        input_output_files: A dict containing input_file: output_file pairs
        recovered_clans_path: A path to the recovered clans project directory
    Returns: None
    """
    for input_file, output_file in input_output_files.items():
        run_clans_headless(recovered_clans_path, input_file, output_file, rounds)


#test
PATH_TO_RECOVERED_CLANS = "/home/aronw/Development/clans-recovered"
input_file = "/home/aronw/Development/Clans-3D/Scores_Evaluation/clans_files/dataset_1_cleaned.clans"
output_file = "/home/aronw/Development/Clans-3D/Scores_Evaluation/clans_files/dataset_1_cleaned_updated.clans"
#run_clans_headless(PATH_TO_RECOVERED_CLANS, input_file, output_file)
