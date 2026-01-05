import sys
import os
sys.path.append(os.path.abspath(".."))
from ToolType import ToolType
from InputFileType import InputFileType
from main import create_clans_file
from recovered_CLANS.utils_old_clans import generate_clans_file_seq_based, run_clans_headless


class ScoresEvaluator:
    def __init__(self, working_dir: str):
        self.working_dir = working_dir
        self.blast_dir = os.path.join(working_dir, "blast_temp")


    def generate_clans_files(self, data: str, input_file_type: InputFileType, tool: ToolType, score: str | None) -> tuple[str, str]:
        if input_file_type == InputFileType.FASTA or input_file_type == InputFileType.TSV or input_file_type == InputFileType.A2M:
            struct_clans_file_path, cleaned_input_file_path = create_clans_file(data, input_file_type, tool, score, structures_dir="structures", out_dir_path=self.working_dir)
            seq_clans_file_path = generate_clans_file_seq_based(cleaned_input_file_path, self.working_dir, self.blast_dir)
            return (struct_clans_file_path, seq_clans_file_path)
        else:
            raise ValueError(f"Invalid input file type: {input_file_type}. Supported types are {InputFileType.FASTA}, {InputFileType.TSV}, and {InputFileType.A2M}.")
        
    
    def cluster_clans_files(self, clans_files: tuple[str, str], path_to_clans_executable: str, rounds_to_cluster: int, p_value: float) -> tuple[str, str]:
        struct_clans_file_clustered = os.path.basename(clans_files[0]).replace(".clans", f"_clustered_r{rounds_to_cluster}_p{p_value}.clans")
        seq_clans_file_clustered = os.path.basename(clans_files[1]).replace(".clans", f"_clustered_r{rounds_to_cluster}_p{p_value}.clans")
        input_output_files = {
            clans_files[0] : f"{self.working_dir}/{struct_clans_file_clustered}",
            clans_files[1] : f"{self.working_dir}/{seq_clans_file_clustered}"
        }
        run_clans_headless(input_output_files, InputFileType.CLANS, path_to_clans_executable, rounds_to_cluster, p_value)
        return (struct_clans_file_clustered, seq_clans_file_clustered)


    def evaluate_clustered_clans_files(self, clustered_clans_files: tuple[str, str]):
        pass