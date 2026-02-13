from clans3d.core.clans_file_generator import ClansFileGenerator
import pandas as pd
from scipy.spatial.distance import pdist
import numpy as np
from clans3d.evaluation.data_normalizer import DataNormalizer


class ClansDataExtractor:
    """Handles extraction and processing of data from CLANS files."""
    
    def __init__(self):
        self.normalizer = DataNormalizer()
    
    def extract_data_from_clans_file_to_df(self, clans_file_path: str) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Extracts scores and coordinates from a CLANS file into pandas DataFrames.

        Args:
            clans_file_path (str): Path to the CLANS file.

        Returns:
            tuple[pd.DataFrame, pd.DataFrame]: A tuple containing:
                - DataFrame of scores with columns ["Sequence_ID_1", "Sequence_ID_2", "Score", "Score_-log10"]
                - DataFrame of coordinates with columns ["Sequence_ID", "x", "y", "z"]
        """
        clans_file_generator = ClansFileGenerator()
        clans_file = clans_file_generator.parse_clans_file(clans_file_path)
        
        df_scores = clans_file.scores.copy()
        df_scores.columns = ["Sequence_ID_1", "Sequence_ID_2", "Score"]
        df_scores["Score_-log10"] = self.normalizer.normalize(df_scores, "Score", "-log10")
        
        df_coord = pd.DataFrame(clans_file.coordinates, columns=["Sequence_ID", "x", "y", "z"])
        
        return df_scores, df_coord
    
    def get_euclidean_from_coordinates(self, df_coordinates: pd.DataFrame) -> pd.DataFrame:
        """
        Computes pairwise Euclidean distances (upper triangle only) from 3D coordinates.
        Also adds another column to the returned DataFrame containing the min-max scaled euclidean distances.

        Args:
            df_coordinates (pd.DataFrame): columns = ["Sequence_ID", "x", "y", "z"]

        Returns:
            pd.DataFrame: columns = ["Sequence_ID_1", "Sequence_ID_2", "euclidean_dist", "euclidean_dist_min_max"]
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
            "euclidean_dist": distances
        })
        
        df_euclidean_dist["euclidean_dist_min_max"] = self.normalizer.normalize(
            df_euclidean_dist, "euclidean_dist", "min-max"
        )
        
        return df_euclidean_dist
