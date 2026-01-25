import sys
import os
sys.path.append(os.path.abspath(".."))
from ToolType import ToolType
from InputFileType import InputFileType
from main import create_clans_file
from recovered_CLANS.utils_old_clans import generate_clans_file_seq_based, run_clans_headless
from ConfigFile import ConfigFile
from ClansFileGenerator import ClansFileGenerator
import pandas as pd
from scipy.spatial.distance import pdist
import numpy as np


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
                            cluster2d: tuple[bool, bool],
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
            "cluster2d": "T" if cluster2d[0] else "F",
            "verbose": int(verbose)
        })
        run_clans_headless(conf_file_struct, path_to_clans_executable)
        
        conf_file_seq = ConfigFile(os.path.join(self.working_dir, seq_clans_file_basename.replace(".clans", ".conf")))
        conf_file_seq.write_config({
            "nographics": "T",
            "load": clans_files[1],
            "dorounds": rounds_to_cluster[1],
            "saveto": seq_clans_file_clustered_path,
            "pval": p_values[1],
            "cluster2d": "T" if cluster2d[1] else "F",
            "verbose": int(verbose)
        })
        run_clans_headless(conf_file_seq, path_to_clans_executable)
        
        return (struct_clans_file_clustered_path, seq_clans_file_clustered_path)


    def evaluate_clustered_clans_files(self, clustered_clans_files: tuple[str, str]):
        print(f"\nEvaluating clustered clans files: {clustered_clans_files[0]} and {clustered_clans_files[1]}")
        
        df_struct_scores, df_struct_coord = self.extract_data_from_clans_file_to_df(clustered_clans_files[0])
        df_struct_euclidean_dist = self.get_pairwise_distances_from_coordinates(df_struct_coord)
        
        df_seq_scores, df_seq_coord = self.extract_data_from_clans_file_to_df(clustered_clans_files[1])
        df_seq_euclidean_dist = self.get_pairwise_distances_from_coordinates(df_seq_coord)
        
        return None


    def extract_data_from_clans_file_to_df(self, clans_file_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        clans_file_generator = ClansFileGenerator()
        clans_file = clans_file_generator.parse_clans_file(clans_file_path)
        df_scores = clans_file.scores
        df_scores.columns = ["Sequence_ID_1", "Sequence_ID_2", "Score"]
        coordinates = clans_file.coordinates
        df_coord = pd.DataFrame(coordinates, columns=["Sequence_ID", "x", "y", "z"])
        return df_scores, df_coord
    
    
    def get_pairwise_distances_from_coordinates(self, df_coordinates: pd.DataFrame) -> pd.DataFrame:
        """
        Computes pairwise Euclidean distances (upper triangle only) from 3D coordinates.

        Args:
            df_coordinates (pd.DataFrame): columns = ["Sequence_ID", "x", "y", "z"]

        Returns:
            pd.DataFrame: columns = ["Sequence_ID_1", "Sequence_ID_2", "euclidean_dist"]
        """
        if df_coordinates.shape[1] != 4:
            raise ValueError("df_coordinates must have exactly 4 columns (Sequence_ID, x, y, z).")
        structures = df_coordinates["Sequence_ID"].to_numpy()
        coords = df_coordinates[["x", "y", "z"]].to_numpy(dtype=float)
        # Compute condensed distance matrix (upper triangle)
        distances = pdist(coords, metric="euclidean")
        # Indices corresponding to upper triangle
        i_idx, j_idx = np.triu_indices(len(structures), k=1)
        df_euclidean_dist = pd.DataFrame({
            "Sequence_ID_1": structures[i_idx],
            "Sequence_ID_2": structures[j_idx],
            "euclidean_dist": distances})
        print(df_euclidean_dist.head())
        return df_euclidean_dist
    