import sys
import os
sys.path.append(os.path.abspath(".."))
from ToolType import ToolType
from InputFileType import InputFileType
from main import create_clans_file
from recovered_CLANS.utils_old_clans import generate_clans_file_seq_based, run_clans_headless
from ConfigFile import ConfigFile
import pandas as pd
from ClansDataExtractor import ClansDataExtractor
from DataNormalizer import DataNormalizer
from ClusterAnalyzer import ClusterAnalyzer
from ClansVisualizer import ClansVisualizer


class ScoresEvaluator:
    """
    Main orchestrator for CLANS file evaluation workflows.
    Provides high-level workflow methods and direct access to specialized component classes.
    
    Component classes (accessible as public attributes):
    - extractor: ClansDataExtractor - Parse CLANS files and extract coordinates
    - normalizer: DataNormalizer - Normalize data columns
    - clustering: ClusterAnalyzer - Clustering and evaluation methods
    - visualizer: ClansVisualizer - Generate visualizations
    """
    
    def __init__(self, working_dir: str):
        self.working_dir = working_dir
        self.blast_dir = os.path.join(working_dir, "blast_temp")
        
        # Expose component classes as public attributes for direct access
        self.extractor = ClansDataExtractor()
        self.normalizer = DataNormalizer()
        self.clustering = ClusterAnalyzer()
        self.visualizer = ClansVisualizer()


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


    def extract_data_from_clans_files(self, clustered_clans_files: tuple[str, str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        High-level workflow method: Extracts and merges data from structural and sequence CLANS files.
        Returns merged DataFrames with _struct and _seq suffixes.
        
        Returns:
            tuple containing:
            1. df_scores: [Sequence_ID_1, Sequence_ID_2, Score_struct, Score_-log10_struct, Score_seq, Score_-log10_seq]
            2. df_euclidean_dist: [Sequence_ID_1, Sequence_ID_2, euclidean_dist_struct, euclidean_dist_min_max_struct, ...]
            3. df_coord: [Sequence_ID, x_struct, y_struct, z_struct, x_seq, y_seq, z_seq]
        
        For direct access to component methods, use:
        - self.extractor.extract_data_from_clans_file_to_df()
        - self.extractor.get_euclidean_from_coordinates()
        """
        print(f"\nEvaluating clustered clans files: {clustered_clans_files[0]} and {clustered_clans_files[1]}")
        # Extract data from structural CLANS file
        df_struct_scores, df_struct_coord = self.extractor.extract_data_from_clans_file_to_df(clustered_clans_files[0])
        df_struct_euclidean_dist = self.extractor.get_euclidean_from_coordinates(df_struct_coord)
        # Extract data from sequence CLANS file
        df_seq_scores, df_seq_coord = self.extractor.extract_data_from_clans_file_to_df(clustered_clans_files[1])
        df_seq_euclidean_dist = self.extractor.get_euclidean_from_coordinates(df_seq_coord)
        # Merge dataframes with suffixes
        df_scores = pd.merge(df_struct_scores, df_seq_scores, on=["Sequence_ID_1", "Sequence_ID_2"], suffixes=("_struct", "_seq"), how="outer")
        df_euclidean_dist = pd.merge(df_struct_euclidean_dist, df_seq_euclidean_dist, on=["Sequence_ID_1", "Sequence_ID_2"], suffixes=("_struct", "_seq"))
        df_coord = pd.merge(df_struct_coord, df_seq_coord, on="Sequence_ID", suffixes=("_struct", "_seq"), how="outer")
        return df_scores, df_euclidean_dist, df_coord

