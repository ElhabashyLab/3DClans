import pandas as pd
from sklearn.cluster import DBSCAN
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from hdbscan import HDBSCAN
import networkx as nx
import igraph as ig
import leidenalg as lg


class ClusterAnalyzer:
    """Handles clustering operations and evaluation metrics for CLANS data."""
    
    # ==================== Clustering Methods ====================
    
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

    def find_clusters_graph_based(
        self, 
        df_scores: pd.DataFrame, 
        scores_col: str, 
        resolution: float = 1.0,
        threshold: float | None = None
    ) -> pd.DataFrame:
        """
        Find clusters in a graph using Leiden community detection.
        The graph is built from similarity scores.

        Args:
            df_scores (pd.DataFrame): DataFrame containing pairwise similarity scores with columns:
                ["Sequence_ID_1", "Sequence_ID_2", <scores_col>]
            scores_col (str): Name of the column containing the similarity scores
            resolution (float): Controls cluster granularity 
                (small -> few big clusters, big -> many small clusters)
            threshold (float | None): Optional minimum score cutoff for edges

        Returns:
            pd.DataFrame: DataFrame with columns ["Sequence_ID", "cluster_id_<scores_type>_Leiden"]
        """
        G = self._build_similarity_graph(df_scores, scores_col, threshold)
        
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

    def _build_similarity_graph(
        self, 
        df_scores: pd.DataFrame, 
        score_col: str, 
        threshold: float | None = None
    ) -> nx.Graph:
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

    # ==================== Evaluation Methods ====================
    
    @staticmethod
    def _validate_columns(df: pd.DataFrame, col1: str, col2: str) -> None:
        """
        Validates that required columns exist in the DataFrame.
        
        Args:
            df (pd.DataFrame): The DataFrame to validate
            col1 (str): First column name
            col2 (str): Second column name
            
        Raises:
            ValueError: If either column is not found
        """
        if col1 not in df.columns or col2 not in df.columns:
            raise ValueError(f"Column '{col1}' or '{col2}' not found in DataFrame.")
    
    @staticmethod
    def _filter_noise_points(df: pd.DataFrame, col1: str, col2: str) -> pd.DataFrame:
        """
        Filters out noise points (cluster label -1) from both clustering results.
        
        Args:
            df (pd.DataFrame): DataFrame containing cluster labels
            col1 (str): First cluster assignment column name
            col2 (str): Second cluster assignment column name
            
        Returns:
            pd.DataFrame: Filtered DataFrame with only valid cluster assignments
            
        Raises:
            ValueError: If no valid cluster assignments remain
        """
        df_valid = df[
            (df[col1] != -1) & 
            (df[col2] != -1)
        ].copy()
        
        if df_valid.empty:
            raise ValueError("No valid cluster assignments found (all points are noise).")
        
        return df_valid
    
    def compute_clustering_agreement(
        self, 
        df_cluster_labels: pd.DataFrame, 
        cluster_col_1: str, 
        cluster_col_2: str
    ) -> dict:
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
        self._validate_columns(df_cluster_labels, cluster_col_1, cluster_col_2)
        
        # Keep only points assigned to a cluster in both clusterings
        df_valid = self._filter_noise_points(df_cluster_labels, cluster_col_1, cluster_col_2)
        
        # Select only the relevant columns for comparison
        df_valid = df_valid[[cluster_col_1, cluster_col_2]]
        
        labels_1 = df_valid[cluster_col_1].to_numpy()
        labels_2 = df_valid[cluster_col_2].to_numpy()

        return {
            "ARI": adjusted_rand_score(labels_1, labels_2),
            "NMI": normalized_mutual_info_score(labels_1, labels_2)
        }
    
    def compute_Jaccard_overlap(
        self, 
        df_cluster_labels: pd.DataFrame, 
        cluster_col_1: str, 
        cluster_col_2: str,
        remove_zero_jaccard: bool = False
    ) -> pd.DataFrame:
        """
        Computes the Jaccard index between all clusters of 2 clustering results.

        Args:
            df_cluster_labels (pd.DataFrame): DataFrame containing cluster labels with columns:
                - Sequence_ID and at least two cluster assignment columns
            cluster_col_1 (str): name of the first cluster assignment column
            cluster_col_2 (str): name of the second cluster assignment column

        Returns:
            pd.DataFrame: DataFrame containing the Jaccard index between all clusters of the two clusterings
                with columns [Cluster_i, Cluster_j, JaccardIndex].
        """
        if cluster_col_1 not in df_cluster_labels.columns or cluster_col_2 not in df_cluster_labels.columns:
            raise ValueError(
                f"Column '{cluster_col_1}' or '{cluster_col_2}' not found in DataFrame."
            )

        # Filter out noise points
        df_valid = df_cluster_labels.loc[
            (df_cluster_labels[cluster_col_1] != -1)
            & (df_cluster_labels[cluster_col_2] != -1),
            ["Sequence_ID", cluster_col_1, cluster_col_2],
        ]
        if df_valid.empty:
            raise ValueError("No valid cluster assignments found (all points are noise).")
        # Assign sequences to clusters
        clusters_1 = {
            cid: set(group["Sequence_ID"])
            for cid, group in df_valid.groupby(cluster_col_1)
        }
        clusters_2 = {
            cid: set(group["Sequence_ID"])
            for cid, group in df_valid.groupby(cluster_col_2)
        }
        # Compute Jaccard index for all cluster pairs
        results = []
        for c1, set_1 in clusters_1.items():
            for c2, set_2 in clusters_2.items():
                intersection = len(set_1 & set_2)
                union = len(set_1 | set_2)
                jaccard_index = intersection / union if union > 0 else 0.0
                if remove_zero_jaccard and jaccard_index == 0.0:
                    continue
                results.append({
                    cluster_col_1: c1,
                    cluster_col_2: c2,
                    "JaccardIndex": jaccard_index,
                })

        return pd.DataFrame(results)
    
    def compute_overlap_coefficient(
        self,
        df_cluster_labels: pd.DataFrame,
        cluster_col_1: str,
        cluster_col_2: str,
        drop_zero: bool = False
    ) -> pd.DataFrame:
        """
        Computes the overlap coefficient (Szymkiewicz-Simpson coefficient) between all clusters 
        of 2 clustering results. The overlap coefficient is defined as:
        
        overlap(A,B) = |A ∩ B| / min(|A|, |B|)
        
        This metric equals 1 when one cluster is a complete subset of the other, making it 
        ideal for detecting subset relationships between clusters.

        Args:
            df_cluster_labels (pd.DataFrame): DataFrame containing cluster labels with columns:
                - Sequence_ID and at least two cluster assignment columns
            cluster_col_1 (str): name of the first cluster assignment column
            cluster_col_2 (str): name of the second cluster assignment column
            drop_zero (bool): If True, exclude cluster pairs with overlap coefficient 0 
                from the output. Default is False.

        Returns:
            pd.DataFrame: DataFrame containing the overlap coefficient between all clusters 
                of the two clusterings with columns [Cluster_i, Cluster_j, OverlapCoefficient, is_smaller].
                The 'is_smaller' column indicates which cluster has fewer sequences:
                - Column name from cluster_col_1 if the first cluster is smaller
                - Column name from cluster_col_2 if the second cluster is smaller
                - "Equal" if both clusters have the same size
        """
        if cluster_col_1 not in df_cluster_labels.columns or cluster_col_2 not in df_cluster_labels.columns:
            raise ValueError(
                f"Column '{cluster_col_1}' or '{cluster_col_2}' not found in DataFrame."
            )

        # Filter out noise points
        df_valid = df_cluster_labels.loc[
            (df_cluster_labels[cluster_col_1] != -1)
            & (df_cluster_labels[cluster_col_2] != -1),
            ["Sequence_ID", cluster_col_1, cluster_col_2],
        ]
        if df_valid.empty:
            raise ValueError("No valid cluster assignments found (all points are noise).")
        
        # Assign sequences to clusters
        clusters_1 = {
            cid: set(group["Sequence_ID"])
            for cid, group in df_valid.groupby(cluster_col_1)
        }
        clusters_2 = {
            cid: set(group["Sequence_ID"])
            for cid, group in df_valid.groupby(cluster_col_2)
        }
        
        # Compute overlap coefficient for all cluster pairs
        results = []
        for c1, set_1 in clusters_1.items():
            for c2, set_2 in clusters_2.items():
                if c1 == c2:
                    continue
                intersection = len(set_1 & set_2)
                min_size = min(len(set_1), len(set_2))
                overlap_coef = intersection / min_size if min_size > 0 else 0.0
                
                if drop_zero and overlap_coef == 0.0:
                    continue
                
                # Determine which cluster is smaller
                if len(set_1) < len(set_2):
                    is_smaller = cluster_col_1
                elif len(set_2) < len(set_1):
                    is_smaller = cluster_col_2
                else:
                    is_smaller = "Equal"
                
                results.append({
                    cluster_col_1: c1,
                    cluster_col_2: c2,
                    "OverlapCoefficient": overlap_coef,
                    "is_smaller": is_smaller
                })

        return pd.DataFrame(results)
    
    
    def compute_cluster_statistics(
        self,
        df_cluster_labels: pd.DataFrame,
        cluster_col: str,
        df_distances: pd.DataFrame,
        distance_col: str,
        df_scores: pd.DataFrame,
        score_col: str
    ) -> pd.DataFrame:
        """
        Computes mean euclidean distance and mean similarity score for each cluster.
        
        Args:
            df_cluster_labels (pd.DataFrame): DataFrame with columns [Sequence_ID, <cluster_col>]
            cluster_col (str): Name of the cluster assignment column
            df_distances (pd.DataFrame): DataFrame with columns [Sequence_ID_1, Sequence_ID_2, <distance_col>]
            distance_col (str): Name of the euclidean distance column
            df_scores (pd.DataFrame): DataFrame with columns [Sequence_ID_1, Sequence_ID_2, <score_col>]
            score_col (str): Name of the similarity score column
        
        Returns:
            pd.DataFrame: DataFrame with columns [Cluster_ID, num_sequences, mean_distance, mean_score]
        """
        if cluster_col not in df_cluster_labels.columns:
            raise ValueError(f"Column '{cluster_col}' not found in cluster labels DataFrame.")
        
        # Filter out noise points
        df_valid = df_cluster_labels[df_cluster_labels[cluster_col] != -1].copy()
        if df_valid.empty:
            raise ValueError("No valid cluster assignments found (all points are noise).")
        
        results = []
        for cluster_id, group in df_valid.groupby(cluster_col):
            seq_ids = set(group["Sequence_ID"])
            num_sequences = len(seq_ids)
            
            # Extract pairwise distances within this cluster
            cluster_distances = df_distances[
                (df_distances["Sequence_ID_1"].isin(seq_ids)) &
                (df_distances["Sequence_ID_2"].isin(seq_ids))
            ]
            mean_distance = cluster_distances[distance_col].mean() if not cluster_distances.empty else None
            # Extract pairwise scores within this cluster
            cluster_scores = df_scores[
                (df_scores["Sequence_ID_1"].isin(seq_ids)) &
                (df_scores["Sequence_ID_2"].isin(seq_ids))
            ]            
            mean_score = cluster_scores[score_col].mean() if not cluster_scores.empty else None
            
            results.append({
                "Cluster_ID": cluster_id,
                "num_sequences": num_sequences,
                "mean_distance": mean_distance,
                "mean_score": mean_score
            })
        
        return pd.DataFrame(results)
