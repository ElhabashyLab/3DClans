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
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.cluster import DBSCAN
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
import networkx as nx
import pandas as pd
import igraph as ig
import leidenalg as lg
from hdbscan import HDBSCAN



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


    def extract_data_from_clans_files(self, clustered_clans_files: tuple[str, str]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Extracts scores and coordinates from a structural and sequence clans file.
        It returns multiple DataFrames:
        
        1. df_scores with columns [Sequence_ID_1, Sequence_ID_2, Score_struct, Score_-log10_struct Score_seq, Score_-log10_seq]
        2. df_euclidean_dist [Sequence_ID_1, Sequence_ID_2, euclidean_dist_struct, euclidean_dist_min_max_struct, euclidean_dist_seq, euclidean_dist_min_max_seq]
        3. df_coord [Sequence_ID, x_struct, y_struct, z_struct, x_seq, y_seq z_seq]
        
        Args:
            clustered_clans_files (tuple[str, str]): Paths to the structural and sequence clans file.

        Returns:
            tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]: df_scores, df_euclidean_dist, df_coord
        """
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
    
    
    def generate_scatter_plot(self, data_x: pd.Series, data_y: pd.Series, 
                        x_label: str|None = None, y_label: str|None = None, 
                        title: str|None = None, save_path: str|None = None):
        """
        Generates a scatter plot comparing two sets of values.

        Args:
            data_x (pd.Series or list-like): Values for the x-axis.
            data_y (pd.Series or list-like): Values for the y-axis.
            x_label (str, optional): Label for x-axis. Defaults to None.
            y_label (str, optional): Label for y-axis. Defaults to None.
            title (str, optional): Plot title. Defaults to None.
            save_path (str, optional): Path to save the figure. If None, plot is shown. Defaults to None.

        Returns:
            None
        """
        data_x = pd.Series(data_x)
        data_y = pd.Series(data_y)
        if len(data_x) != len(data_y):
            raise ValueError("data_x and data_y must have the same length.")

        # Default labels if not provided
        if x_label is None:
            x_label = "X-axis"
        if y_label is None:
            y_label = "Y-axis"
        if title is None:
            title = f"{y_label} vs {x_label}"

        plt.figure(figsize=(6, 6))
        sns.scatterplot(x=data_x, y=data_y, alpha=0.7, edgecolor=None)
        
        # Plot identity line y=x
        min_val = min(data_x.min(), data_y.min())
        max_val = max(data_x.max(), data_y.max())
        plt.plot([min_val, max_val], [min_val, max_val], linestyle='--', color='gray', alpha=0.5)

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300)
            plt.close()
        else:
            plt.show()
        
        
    def find_clusters_density_based(
        self, 
        df_coord: pd.DataFrame,
        coord_type: str, 
        algorithm: str = "HDBSCAN",
        eps: float = 0.5, 
        min_samples: int = 5,
        min_cluster_size: int = 5
    ) -> pd.DataFrame:
        """
        Performs density-based clustering on 3D CLANS coordinates using DBSCAN or HDBSCAN.

        Args:
            df_coord (pd.DataFrame):
                DataFrame containing coordinates. Must include columns:
                - Sequence_ID, x_<coord_type>, y_<coord_type>, z_<coord_type>
            coord_type (str):
                Identifier for which coordinate set to cluster.
                Example: "struct" or "seq"
            algorithm (str):
                Clustering algorithm to use: "DBSCAN" or "HDBSCAN".
                Default is "HDBSCAN" (recommended for most cases).
            eps (float):
                Maximum distance between two points to be considered neighbors.
                Only used for DBSCAN. Default is 0.5.
            min_samples (int):
                Minimum number of points required to form a dense region (DBSCAN)
                or minimum samples in neighborhood (HDBSCAN).
                Default is 5.
            min_cluster_size (int):
                Minimum number of points required to form a cluster.
                Only used for HDBSCAN. Default is 5.

        Returns:
            pd.DataFrame:
                DataFrame with columns: [Sequence_ID, cluster_id_<coord_type>_<algorithm>]
                Noise points are labeled as -1.

        Raises:
            ValueError: If required columns are missing or invalid algorithm specified.
        """
        required_cols = ["Sequence_ID", f"x_{coord_type}", f"y_{coord_type}", f"z_{coord_type}"]
        missing = [c for c in required_cols if c not in df_coord.columns]
        if missing:
            raise ValueError(
                f"Missing required columns for coord_type='{coord_type}': {missing}"
            )
        coords = df_coord[required_cols[1:]].to_numpy(dtype=float)
        
        if algorithm.upper() == "DBSCAN":
            clusterer = DBSCAN(eps=eps, min_samples=min_samples)
            labels = clusterer.fit_predict(coords)
            algo_name = "DBSCAN"
        elif algorithm.upper() == "HDBSCAN":
            clusterer = HDBSCAN(
                min_cluster_size=min_cluster_size,
                min_samples=min_samples
            )
            labels = clusterer.fit_predict(coords)
            algo_name = "HDBSCAN"
        else:
            raise ValueError(
                f"Invalid algorithm '{algorithm}'. Must be 'DBSCAN' or 'HDBSCAN'."
            )
        df_out = pd.DataFrame({
            "Sequence_ID": df_coord["Sequence_ID"],
            f"cluster_id_{coord_type}_{algo_name}": labels
        })
        return df_out


    def build_similarity_graph(self, df_scores: pd.DataFrame, score_col: str, threshold: float | None = None) -> nx.Graph:
        """
        Builds a weighted similarity graph from pairwise scores.

        Args:
            df_scores (pd.DataFrame): columns = [Sequence_ID_1, Sequence_ID_2, score_col]
            score_col (str): column with similarity score (Bigger should mean more similar)
            threshold (float): optional minimum score cutoff

        Returns:
            networkx.Graph
        """
        G = nx.Graph()
        for _, row in df_scores.iterrows():
            weight = row[score_col]
            if threshold is not None:
                if weight < threshold:
                    continue
            G.add_edge(
                row["Sequence_ID_1"],
                row["Sequence_ID_2"],
                weight=weight
            )
        return G


    def find_clusters_graph_based(self, df_scores: pd.DataFrame, scores_col: str, resolution: float = 1.0) -> pd.DataFrame:
        """
        Find clusters in a graph using Leiden community detection.
        The graph is built from similarity scores.

        Args:
            df_scores (pd.DataFrame): DataFrame containing pairwise similarity scores with columns:
                ["Sequence_ID_1", "Sequence_ID_2", Score_<scores_type>]
            scores_col (str): Name of the column containing the similarity scores
            resolution (float): Controls cluster granularity 
                (small -> few big clusters, big -> many small clusters)

        Returns:
            pd.DataFrame: DataFrame with columns ["Sequence_ID", "cluster_id_<scores_type>"]
        """
        G = self.build_similarity_graph(df_scores, scores_col)
        # Convert NetworkX graph to igraph
        g_ig = ig.Graph.from_networkx(G)
        g_ig.vs["name"] = list(G.nodes)
        # Run Leiden community detection
        partition = lg.find_partition(
            g_ig,
            lg.RBConfigurationVertexPartition,
            weights="weight",
            resolution_parameter=resolution
        )
        clusters_dict = {g_ig.vs[i]["name"]: cid for i, cid in enumerate(partition.membership)}
        cluster_id_col_name = f"cluster_id_{scores_col}_Leiden"
        df_out = pd.DataFrame.from_dict(clusters_dict, orient="index", columns=[cluster_id_col_name])
        df_out.index.name = "Sequence_ID"
        df_out.reset_index(inplace=True)
        return df_out


    def compute_clustering_agreement(self, df_cluster_labels: pd.DataFrame, 
                                     cluster_col_1: str, cluster_col_2: str) -> dict:
        """
        Computes clustering agreement metrics (ARI and NMI) between two clustering results.
        
        Args:
            df_cluster_labels (pd.DataFrame): DataFrame containing cluster labels with columns:
                - Sequence_ID and at least two cluster assignment columns
            cluster_col_1 (str): Name of the first cluster assignment column
            cluster_col_2 (str): Name of the second cluster assignment column
        
        Returns:
            dict: Dictionary containing:
                - 'ARI': Adjusted Rand Index (range: -1 to 1, higher is better)
                - 'NMI': Normalized Mutual Information (range: 0 to 1, higher is better)
        
        Raises:
            ValueError: If specified columns are not found in the DataFrame
        """
        if cluster_col_1 not in df_cluster_labels.columns or cluster_col_2 not in df_cluster_labels.columns:
            raise ValueError(f"Column '{cluster_col_1}' or '{cluster_col_2}' not found in DataFrame.")
        
        # Keep only points assigned to a cluster in both clusterings
        df_valid = df_cluster_labels.loc[
            (df_cluster_labels[cluster_col_1] != -1) & (df_cluster_labels[cluster_col_2] != -1),
            [cluster_col_1, cluster_col_2]]

        if df_valid.empty:
            raise ValueError("No valid cluster assignments found (all points are noise).")

        labels_1 = df_valid[cluster_col_1].to_numpy()
        labels_2 = df_valid[cluster_col_2].to_numpy()

        return {
            "ARI": adjusted_rand_score(labels_1, labels_2),
            "NMI": normalized_mutual_info_score(labels_1, labels_2)
        }
    