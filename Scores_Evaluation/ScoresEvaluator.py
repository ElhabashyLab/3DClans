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


    def extract_data_from_clans_files(self, clustered_clans_files: tuple[str, str]):
        print(f"\nEvaluating clustered clans files: {clustered_clans_files[0]} and {clustered_clans_files[1]}")
        # get data from struct clans file
        df_struct_scores, df_struct_coord = self.extract_data_from_clans_file_to_df(clustered_clans_files[0])
        df_struct_euclidean_dist = self.get_euclidean_from_coordinates(df_struct_coord)
        # get data from seq clans file
        df_seq_scores, df_seq_coord = self.extract_data_from_clans_file_to_df(clustered_clans_files[1])
        df_seq_euclidean_dist = self.get_euclidean_from_coordinates(df_seq_coord)
        # combine dataframes
        df_scores = pd.merge(df_struct_scores, df_seq_scores, on=["Sequence_ID_1", "Sequence_ID_2"], suffixes=("_struct", "_seq"), how="outer")
        df_euclidean_dist = pd.merge(df_struct_euclidean_dist, df_seq_euclidean_dist, on=["Sequence_ID_1", "Sequence_ID_2"], suffixes=("_struct", "_seq"))
        df_coord = pd.merge(df_struct_coord, df_seq_coord, on="Sequence_ID", suffixes=("_struct", "_seq"), how="outer")
        return df_scores, df_euclidean_dist, df_coord


    def extract_data_from_clans_file_to_df(self, clans_file_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Extracts scores and coordinates from a CLANS file into pandas DataFrames.

        Args:
            clans_file_path (str): Path to the CLANS file.

        Returns:
            tuple[pd.DataFrame, pd.DataFrame]: A tuple containing:
                - DataFrame of scores with columns ["Sequence_ID_1", "Sequence_ID_2", "Score"]
                - DataFrame of coordinates with columns ["Sequence_ID", "x", "y", "z"]
        """
        clans_file_generator = ClansFileGenerator()
        clans_file = clans_file_generator.parse_clans_file(clans_file_path)
        df_scores = clans_file.scores
        df_scores.columns = ["Sequence_ID_1", "Sequence_ID_2", "Score"]
        df_scores["Score_-log10"] = self.normalize_df_column(df_scores, "Score", "-log10")
        coordinates = clans_file.coordinates
        df_coord = pd.DataFrame(coordinates, columns=["Sequence_ID", "x", "y", "z"])
        return df_scores, df_coord
    
    
    def get_euclidean_from_coordinates(self, df_coordinates: pd.DataFrame) -> pd.DataFrame:
        """
        Computes pairwise Euclidean distances (upper triangle only) from 3D coordinates.
        Also adds another column to the returned DataFrame containing the min-max scaled euclidean distances.

        Args:
            df_coordinates (pd.DataFrame): columns = ["Sequence_ID", "x", "y", "z"]

        Returns:
            pd.DataFrame: columns = ["Sequence_ID_1", "Sequence_ID_2", "euclidean_dist, euclidean_dist_min_max"]
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
        df_euclidean_dist["euclidean_dist_min_max"] = self.normalize_df_column(df_euclidean_dist, "euclidean_dist", "min-max") 
        return df_euclidean_dist
    
    
    def normalize_df_column(self, df: pd.DataFrame, column_name: str, normalization_type: str) -> pd.Series:
        """
        Normalizes a specified column in the DataFrame using the given normalization type.

        Args:
            df (pd.DataFrame): The DataFrame containing the column to normalize.
            column_name (str): The name of the column to normalize.
            normalization_type (str): The type of normalization to apply. Supported types: ("min-max", "z-score", "-log10").

        Returns:
            pd.Series: A Series with the normalized column.
        """
        if column_name not in df.columns:
            raise ValueError(f"Column '{column_name}' not found in DataFrame.")
        df = df.copy()
        col = df[column_name].astype(float)
        normalizers = {
            "min-max": self._normalize_min_max,
            "z-score": self._normalize_z_score,
            "-log10": self._normalize_log10
        }
        if normalization_type not in normalizers:
            raise ValueError(
                f"Unsupported normalization_type '{normalization_type}'. "
                "Supported: 'min-max', 'z-score', '-log10'."
            )
        series = normalizers[normalization_type](col)
        return series
    
    
    def _normalize_min_max(self, col: pd.Series) -> pd.Series:
        """
        Applies min-max normalization to a pandas Series.
        Scales all values to the range [0, 1]. If the column is constant
        (max == min), all values are set to 0.

        Args:
            col (pd.Series): The input numeric column to normalize.

        Returns:
            pd.Series: A Series with min-max normalized values in [0, 1].
        """
        min_val = col.min()
        max_val = col.max()
        if max_val == min_val:
            return pd.Series(0.0, index=col.index)
        return (col - min_val) / (max_val - min_val)


    def _normalize_z_score(self, col: pd.Series) -> pd.Series:
        """
        Applies z-score normalization (standardization) to a pandas Series.
        Each value is transformed as (x - mean) / std. If the column has zero
        variance (std == 0), all values are set to 0.

        Args:
            col (pd.Series): The input numeric column to normalize.

        Returns:
            pd.Series: A Series with z-score normalized values.
        """
        std = col.std()
        if std == 0:
            return pd.Series(0.0, index=col.index)
        return (col - col.mean()) / std


    def _normalize_log10(self, col: pd.Series) -> pd.Series:
        """
        Applies -log10 transformation to a pandas Series making small values more interpretable.
        Raises a ValueError if the column contains zero or negative values.

        Args:
            col (pd.Series): The input numeric column to transform. Must be > 0.

        Returns:
            pd.Series: A Series where each value is transformed as -log10(x).

        Raises:
            ValueError: If any value in the column is <= 0.
        """
        if (col <= 0).any():
            raise ValueError("Cannot apply -log10 to zero or negative values.")
        return pd.Series(-np.log10(col.values), index=col.index)