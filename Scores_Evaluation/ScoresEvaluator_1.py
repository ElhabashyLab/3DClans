import sys
import os
sys.path.append(os.path.abspath(".."))
from ToolType import ToolType
from InputFileType import InputFileType
from main import create_clans_file
from recovered_CLANS.utils_old_clans import generate_clans_file_seq_based, run_clans_headless
from ConfigFile import ConfigFile


class ScoresEvaluator:
    def __init__(self, working_dir: str):
        self.working_dir = working_dir
        self.blast_dir = os.path.join(working_dir, "blast_temp")


    def generate_clans_files(self, data: str, input_file_type: InputFileType, tool: ToolType, score: str | None) -> tuple[str, str]:
        if input_file_type == InputFileType.FASTA or input_file_type == InputFileType.TSV or input_file_type == InputFileType.A2M:
            struct_clans_file_path, cleaned_input_file_as_fasta_path = create_clans_file(data, input_file_type, tool, score, structures_dir="structures", out_dir_path=self.working_dir)
            seq_clans_file_path = generate_clans_file_seq_based(cleaned_input_file_as_fasta_path, self.working_dir, self.blast_dir)
            return (struct_clans_file_path, seq_clans_file_path)
        else:
            raise ValueError(f"Invalid input file type: {input_file_type}. Supported types are {InputFileType.FASTA}, {InputFileType.TSV}, and {InputFileType.A2M}.")
        
    
    def cluster_clans_files(self,
                            path_to_clans_executable: str,
                            clans_files: tuple[str, str],
                            rounds_to_cluster: tuple[int, int],
                            p_values: tuple[float, float],
                            cluster2D: tuple[bool, bool],
                            verbose: bool) -> tuple[str, str]:
        struct_clans_file_basename = os.path.basename(clans_files[0])
        seq_clans_file_basename = os.path.basename(clans_files[1])
        struct_clans_file_clustered_name = struct_clans_file_basename.replace(".clans", f"_clustered_r_{rounds_to_cluster[0]}_p_{p_values[0]}.clans")
        seq_clans_file_clustered_name = seq_clans_file_basename.replace(".clans", f"_clustered_r_{rounds_to_cluster[1]}_p_{p_values[1]}.clans")
        struct_clans_file_clustered_path = os.path.join(self.working_dir, struct_clans_file_clustered_name)
        seq_clans_file_clustered_path =  os.path.join(self.working_dir, seq_clans_file_clustered_name)

        conf_file_struct = ConfigFile(os.path.join(self.working_dir, struct_clans_file_basename.replace(".clans", ".conf")))
        conf_file_struct.write_config({
            "nographics": "T",
            "load": clans_files[0],
            "dorounds": rounds_to_cluster[0],
            "saveto": struct_clans_file_clustered_path,
            "pval": p_values[0],
            "verbose": int(verbose)
        })
        conf_file_seq = ConfigFile(os.path.join(self.working_dir, seq_clans_file_basename.replace(".clans", ".conf")))
        conf_file_seq.write_config({
            "nographics": "T",
            "load": clans_files[1],
            "dorounds": rounds_to_cluster[1],
            "saveto": seq_clans_file_clustered_path,
            "pval": p_values[1],
            "verbose": int(verbose)
        })

        run_clans_headless(
            conf_file_struct,
            cluster2D[0],
            path_to_clans_executable)

        run_clans_headless(
            conf_file_seq,
            cluster2D[1],
            path_to_clans_executable)
        
        return (struct_clans_file_clustered_path, seq_clans_file_clustered_path)


    def evaluate_clustered_clans_files(self, clustered_clans_files: tuple[str, str]):
        pass