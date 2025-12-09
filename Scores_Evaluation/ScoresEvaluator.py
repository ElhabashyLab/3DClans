import sys
import os
sys.path.append(os.path.abspath(".."))
from Dataset_Generator.DatasetGenerator import DatasetGenerator
from old_clans.utils_old_clans import run_clans_headless
from utils_for_PDB import copy_dir_content, reset_dir_content, generate_fasta_from_uids_with_regions, fetch_pdbs
from StructSimComputer import StructSimComputer
from ClansFileGenerator import ClansFileGenerator
from ClansFile import ClansFile
from ToolType import ToolType
from InputFileType import InputFileType
import pandas as pd 
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import networkx as nx
from typing import Optional
from skbio.stats.distance import mantel
from skbio import DistanceMatrix
from sklearn.cluster import DBSCAN
from sklearn.metrics import adjusted_rand_score, normalized_mutual_info_score
from scipy.spatial import procrustes


class ScoresEvaluator:
    """
    This class is responsible for evaluating the structure similarity scores by comparing them to sequence similarity scores.
    The datasets used for the evaluation can be generated with the DatasetGenerator class with the given seeds.
    Otherwise the user can provide existing datasets.
    With the data, the structure similarity scores and sequence similarity scores are computed and clans files are generated for each.
    The comparison is done with the scores as well as with inferred graphs and clusters.
    
    Order of method calls for evaluation:
    1. __init__: Initialize the ScoresEvaluator with necessary parameters.
    2. _initialize_evaluation: provide or generate the data, compute scores, and generate clans files.
    3. evaluate: Perform the evaluation on the generated clans files and return the results.
    4. _display_evaluation_results: Display the evaluation results including plots and dataframes
    """
    def __init__(self, path_to_recovered_clans: str):
        """
        Innitializes the ScoresEvaluator by setting up necessary directories and creating instances of required classes.

        Args:
        path_to_recovered_clans (str): Path to the recovered clans.jar file.
        rounds_to_cluster (int): Number of rounds to use for clustering in recovered clans.
        generate_datasets (bool): Whether to generate datasets or use existing data.
        data_for_datasets (list): A list containing number_of_datasets (int), size_of_datasets (int) and seeds (list of str) for dataset generation.
        """
        self.working_dir = os.path.dirname(os.path.abspath("ScoresEvaluator.py"))
        self._set_up_dirs()
        self.dataset_generator = DatasetGenerator(self.datasets_dir)
        self.scores_computer = StructSimComputer(foldseek_score="evalue")
        self.tool_type = ToolType.FOLDSEEK
        self.struct_clans_generator = ClansFileGenerator(self.clans_files_structsim_dir)
        self.seq_clans_generator = ClansFileGenerator(self.clans_files_seqsim_dir) 
        self.path_to_recovered_clans = path_to_recovered_clans


    def initialize_evaluation(self,
                              rounds_to_cluster: int,
                              datasets_meta_data: Optional[dict] = None,
                              use_existing_dataset: bool = False,
                              path_to_dir_of_existing_datasets: Optional[str] = None,
                              datasets_file_type: InputFileType = InputFileType.FASTA
                              ) -> dict:
        """
        Generates or uses existing datasets.
        Downloads the corresponding PDB files of the datasets, computes structure similarity scores and creates clans files from those.
        After, it runs the old recovered clans 'rounds_to_cluster'-times on the generated clans files and saves the updated clans files.
        This function should be called before running the actual evaluation.
        Args:
            rounds_to_cluster: Number of rounds to use for clustering in recovered clans.
            datasets_meta_data: A dictionary containing for each dataset; size_of_dataset (int), number of clusters (int) and seeds (list of str).
            use_existing_dataset: Whether to use existing datasets or generate new ones.
            path_to_dir_of_existing_datasets: Path to existing datasets if use_existing_dataset is True.
            datasets_file_type: The type of the dataset files in the existing datasets. If None, FASTA is assumed.
        Retruns:
            dict: A dictionary mapping structural clans file paths to sequence clans file paths.
        """
        print("Initializing evaluation...")
        self._set_up_datasets_dir(use_existing_dataset, path_to_dir_of_existing_datasets, datasets_meta_data)
        uids_with_regions_for_each_dataset = self._download_structures(self.datasets_dir, self.structures_dir, datasets_file_type)
        paths_to_cleaned_datasets = self._generate_fasta_from_uids_with_regions_for_each_dataset(uids_with_regions_for_each_dataset, self.datasets_dir, datasets_file_type)
        scores_for_each_dataset = self._compute_struct_scores(self.structures_dir)
        self._compute_struct_clans_files(scores_for_each_dataset, paths_to_cleaned_datasets)
        
        print("Running recovered clans.jar on the generated structural clans files...")
        input_output_dict_structural = self._generate_input_output_files_dict(self.clans_files_structsim_dir, self.clans_files_structsim_dir)
        run_clans_headless(self.path_to_recovered_clans, input_output_dict_structural, input_file_type=InputFileType.CLANS, rounds=rounds_to_cluster)
        
        print("Running recovered clans.jar with fasta files (using sequence similarity)...")
        input_output_dict_sequence = self._generate_input_output_files_dict(self.datasets_dir, self.clans_files_seqsim_dir)
        run_clans_headless(self.path_to_recovered_clans, input_output_dict_sequence, input_file_type=InputFileType.FASTA, blast_dir=self.blast_dir, clans_generator=self.seq_clans_generator, rounds=rounds_to_cluster)
        print("Evaluation initialized. Clans files generated and clustered for each dataset with structure similarity scores and sequence similarity scores.")
        return self._match_clans_files_for_comparison(self.clans_files_structsim_dir, self.clans_files_seqsim_dir)
    
    
    def evaluate(self, structural_to_sequence_clans_files: dict):
        """
        Args:
            structural_to_sequence_clans_files: A dictionary mapping structural clans file paths to sequence clans file paths.
        """
        evaluation_results = []
        figures = []
        
        for struct_clans_file, seq_clans_file in structural_to_sequence_clans_files.items():
            
            print(f"Evaluating structural clans file {struct_clans_file} with sequence clans file {seq_clans_file}...")
            parsed_struct_clans_file = ClansFileGenerator.parse_clans_file(struct_clans_file)
            parsed_seq_clans_file = ClansFileGenerator.parse_clans_file(seq_clans_file)
            
            # computing metrics
            scores_combined_df, coordinates_combined_df = self._prepare_scores_and_coordinates(parsed_struct_clans_file, parsed_seq_clans_file)
            print(scores_combined_df)
            print(coordinates_combined_df)
            
            # numerical evaluation
            figures.append(self._plot_similarity_scatter(scores_combined_df, x_col="score_struct_log10_min-max-scaled", y_col="score_seq_log10_min-max-scaled", prevent_display=False))
            df_numerical_comparison = self._compare_numerically(scores_combined_df)
            print(df_numerical_comparison)
            
            # overall graph evaluation
            coordinates_struct = coordinates_combined_df[["PDBchain1", "x_struct", "y_struct", "z_struct"]]
            coordinates_seq = coordinates_combined_df[["PDBchain1", "x_seq", "y_seq", "z_seq"]]
            G_struct = self._build_graph_from_scores(scores_combined_df[["PDBchain1", "PDBchain2", "score_struct"]], coordinates_struct)
            G_seq = self._build_graph_from_scores(scores_combined_df[["PDBchain1", "PDBchain2", "score_seq"]], coordinates_seq)
            df_graph_comparison = self._compare_graphs(G_struct, G_seq)
            print(df_graph_comparison)
            
            # cluster evaluation
            coord_and_cluster_labels_combined_df = self._get_cluster_labels(coordinates_combined_df)
            cluster_labels_df = coord_and_cluster_labels_combined_df[["PDBchain1", "cluster_label_struct", "cluster_label_seq"]]
            ari = adjusted_rand_score(cluster_labels_df["cluster_label_struct"], cluster_labels_df["cluster_label_seq"])
            nmi = normalized_mutual_info_score(cluster_labels_df["cluster_label_struct"], cluster_labels_df["cluster_label_seq"])
            coord_struct_z_scored = coord_and_cluster_labels_combined_df[["x_struct_z-scores", "y_struct_z-scores", "z_struct_z-scores"]]
            coord_seq_z_scored = coord_and_cluster_labels_combined_df[["x_seq_z-scores", "y_seq_z-scores", "z_seq_z-scores"]]
            mtx1, mtx2, disparity = procrustes(coord_seq_z_scored, coord_struct_z_scored)
            df_cluster_comparison = pd.DataFrame({"ari": [ari], "nmi": [nmi], "procrustes_disparity": [disparity]})
            df_metrics_per_cluster, df_cluster_overlap = self._compare_clusters(scores_combined_df, cluster_labels_df)
            print(cluster_labels_df)
            print(df_cluster_comparison)
            print(df_metrics_per_cluster)
            print(df_cluster_overlap)
        
        print("Evaluation completed.")
        return evaluation_results, figures
        
        
    def display_evaluation_results(self, results_per_dataset, figures):
        print("Scatterplots: ")
        for fig in figures:
            fig
        print("Numerical comparison results:")
        print(results_per_dataset)
        # add further results here:    
    
    
    def _set_up_datasets_dir(self, use_existing_dataset, path_to_dir_of_existing_datasets, datasets_meta_data):
        """
        Generates datasets with the given datasets_meta_data in self.datasets_dir.
        If use_existing_dataset is True, it copys the content of the folder of existing datasets to self.datasets_dir.

        Args:
            use_existing_dataset (bool): Whether to use an existing dataset directory.
            path_to_dir_of_existing_datasets (str): Path to existing datasets if use_existing_dataset is True.
            datasets_meta_data (dict): Metadata for generating datasets if not using existing datasets.

        Raises:
            ValueError: If use_existing_dataset is True but no path is provided.
        """
        if use_existing_dataset:
            if path_to_dir_of_existing_datasets is None:
                raise ValueError("Path to existing datasets must be provided if use_existing_dataset is True.")
            copy_dir_content(path_to_dir_of_existing_datasets, self.datasets_dir)
        else:
            if datasets_meta_data is None:
                raise ValueError("Datasets metadata must be provided if not using existing datasets.")
            self._generate_datasets(datasets_meta_data)
            
        
    def _set_up_dirs(self, leave_as_is: list = []):
        """
        Sets up the necessary directories for the evaluation.
        Creates directories for datasets, structures, clans files and blast files.
        Deletes the content of the directories if they already exist.
        Args:
            leave_as_is (list): A list of directory paths which should not be deleted if they already exist
        """
        self.datasets_dir: str = f"{self.working_dir}/datasets"
        self.structures_dir = f"{self.working_dir}/structures"
        self.clans_files_structsim_dir = f"{self.working_dir}/clans_files_structsim" # dir for clans files generated with structure similarity scores
        self.clans_files_seqsim_dir = f"{self.working_dir}/clans_files_seqsim" # dir for clans files generated with sequence similarity scores
        self.blast_dir = f"{self.working_dir}/blast_files"
        dirs = [self.datasets_dir, self.structures_dir, self.clans_files_structsim_dir, self.blast_dir, self.clans_files_seqsim_dir]
        for dir_path in dirs:
            if dir_path not in leave_as_is:
                reset_dir_content(dir_path)


    def _generate_datasets(self, dataset_information: dict):
        """
        Generates datasets (fasta_files) by calling the DatasetGenerator class.
        The datasets are saved in the out_dir.
        
        Args:
            dataset_information: A dictionary containing for each dataset; size_of_dataset (int), number of clusters (int) and seeds (list of str).
            out_dir: directory where to save the datasets
        """
        for dataset_name, meta_data in dataset_information.items():
            size = meta_data["size_of_dataset"]
            n_cl = meta_data["number_of_clusters"]
            seeds = meta_data["seeds"]
            print(f"Generating dataset {dataset_name} with {size} sequences and {n_cl} clusters...")
            self.dataset_generator.generate(size, n_cl, seeds, f"{dataset_name}.fasta")
            
    
    def _prepare_scores_and_coordinates(self, struct_clans_file: ClansFile, seq_clans_file: ClansFile) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Prepares the scores and coordinates dataframes for comparison by log-transforming and normalizing the scores and normalizing the coordinates.

        Args:
            struct_clans_file (ClansFile): A ClansFile Object based on structural similarity.
            seq_clans_file (ClansFile): A ClansFIle Object based on sequence similarity.

        Returns:
            tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two DataFrames: the combined scores DataFrame and the combined coordinates DataFrame.
        """
        scores_combined_df, coordinates_combined_df = self._get_combined_scores_and_coordinates_df(struct_clans_file, seq_clans_file)
        scores_combined_df = self._log_transform(scores_combined_df, ["score_struct", "score_seq"])
        scores_combined_df = self._normalize_scores(scores_combined_df, ["score_struct_log10", "score_seq_log10"], "z-score")
        scores_combined_df = self._normalize_scores(scores_combined_df, ["score_struct_log10", "score_seq_log10"], "min-max")
        coordinates_combined_df = self._normalize_scores(coordinates_combined_df, ["x_struct", "y_struct", "z_struct", "x_seq", "y_seq", "z_seq"], "z-score")
        return scores_combined_df, coordinates_combined_df
        
        
    def _normalize_scores(self, df: pd.DataFrame, columns: list, norm_type: str) -> pd.DataFrame: 
        """
        Normalizes the given columns of the given dataframe with the specified normalization-types.

        Args:
            df: Dataframe containing columns with scores.
            columns: Names of columns which should be normalized.
            norm_type: Type of normalization for each column in columns. Possible values are "min-max" and "z-score"

        Returns:
            The df with normalized scores.
        """
        df_norm = df.copy()
        for col in columns:
            if col not in df_norm.columns:
                raise ValueError(f"Column '{col}' not found in DataFrame.")
            series = df_norm[col]
            if norm_type == "min-max":
                df_norm[col + "_min-max-scaled"] = self._normalize_min_max(series)
            elif norm_type == "z-score":
                df_norm[col + "_z-scores"] = self._normalize_to_z_scores(series)
            else:
                raise ValueError(f"Unsupported normalization type '{norm_type}'. "
                                "Use 'min-max' or 'z-score'.")
        return df_norm


    def _log_transform(self, df: pd.DataFrame, columns: list, base: float = 10) -> pd.DataFrame:
        """
        Log-transforms the specified columns of a Pandas DataFrame using a given logarithm base.
        Automatically handles zeros and non-positive values by replacing them with the smallest positive float.

        Args:
            df (pd.DataFrame): The DataFrame containing columns to be transformed.
            columns (list): List of column names to be log-transformed.
            base (float): The logarithm base (default: 10).

        Returns:
            pd.DataFrame: A copy of the input DataFrame with additional log-transformed columns.
                        Each transformed column is named '<original>_log<base>'.
        """
        df_log = df.copy()
        for col in columns:
            if col not in df_log.columns:
                raise ValueError(f"Column '{col}' not found in DataFrame.")
            # Replace zeros or negatives with smallest positive float to avoid log(0)
            series = df_log[col]
            safe_series = series.copy()
            safe_series[safe_series <= 0] = np.nextafter(0, 1)
            if base == 10:
                transformed = -np.log10(safe_series)
            elif base == np.e:
                transformed = -np.log(safe_series)
            else:
                transformed = -np.log(safe_series) / np.log(base)

            # Add new column
            new_col_name = f"{col}_log{base}"
            df_log[new_col_name] = pd.Series(transformed, index=df_log.index, name=new_col_name)
        return df_log
    
    
    def _normalize_min_max(self, scores: pd.Series) -> pd.Series:
        """
        Normalizes a given Pandas Series with min-max scaling.

        Args:
            scores (pd.Series): A Series containing scores.

        Returns:
            scores_norm: A Series containing normalized (min-max scaled) scores.
        """
        min_val, max_val = scores.min(), scores.max()
        if pd.isna(min_val) or pd.isna(max_val) or min_val == max_val:
            scores_norm = scores # avoid divide-by-zero
        else:
            scores_norm = (scores - min_val) / (max_val - min_val)
        return scores_norm
        
        
    def _normalize_to_z_scores(self, scores: pd.Series) -> pd.Series:
        """
        Normalizes a given Pandas Series to z-scores.

        Args:
            scores (pd.Series): A Series containing scores.

        Returns:
            scores_norm: A Series containing normalized (z-score) scores.
        """
        mean, std = scores.mean(), scores.std()
        if pd.isna(std) or std == 0:
            scores_norm = scores
        else:
            scores_norm = (scores - mean) / std
        return scores_norm
    
    
    def _plot_similarity_scatter(self, df, x_col='score_struct', y_col='score_seq', title='Structural vs Sequence Similarity', prevent_display=False):
        """
        Create a scatter plot with a regression line between two columns of a DataFrame.
        Args:
            df (pd.DataFrame): The DataFrame containing the data.
            x_col (str): The name of the column to be plotted on the x-axis.
            y_col (str): The name of the column to be plotted on the y-axis.
            title (str): The title of the plot.
        Returns:
            matplotlib.figure.Figure: The generated scatter plot with regression line which can be displayed or saved.
        """
        fig, ax = plt.subplots(figsize=(6, 6))
        sns.regplot(
            data=df,
            x=x_col,
            y=y_col,
            scatter_kws={'alpha': 0.6, 's': 40},
            line_kws={'color': 'red'},
            ax=ax
        )
        ax.set_title(title)
        ax.set_xlabel(x_col.replace('_', ' ').title())
        ax.set_ylabel(y_col.replace('_', ' ').title())
        ax.grid(True)
        if prevent_display:
            plt.close(fig)  # prevent immediate display in notebooks
        return fig
        

    def _compare_numerically(self, df_combined: pd.DataFrame) -> pd.DataFrame:
        """
        Compares the scores of columns contained in df_combined numerically.
        Returns a dataframe with the comparison results.
        The following metrics are computed: 
        Correlation(Spearman and Pearson), RMSD, mantel_corr and mantel_p_value
        Args:
            df_combined: A dataframe with different structural and sequence scores.
        Returns:
            df_results: A dataframe containg Correlation(Spearman and Pearson), RMSD values, mantel_corr and mantel_p_value.
        """
        df_combined = df_combined.apply(pd.to_numeric)
        spearman_corr = df_combined["score_struct_log10"].corr(df_combined["score_seq_log10"], method='spearman')
        pearson_corr = df_combined["score_struct_log10_z-scores"].corr(df_combined["score_seq_log10_z-scores"], method='pearson')
        # rmsd requires similarity scores and not distance scores -> 1 - distance
        df_combined["score_struct_log10_sim"] = -df_combined["score_struct_log10"]
        df_combined["score_seq_log10_sim"] = -df_combined["score_seq_log10"]
        rmsd = ((df_combined["score_struct_log10_sim"] - df_combined["score_seq_log10_sim"]) ** 2).mean() ** 0.5
        dist_mat_struct = self._convert_df_to_distance_matrix(df_combined, "score_struct_log10")
        dist_mat_seq = self._convert_df_to_distance_matrix(df_combined, "score_seq_log10")
        mantel_corr, mantel_p_value, _ = mantel(dist_mat_struct, dist_mat_seq, method='spearman', permutations=999)
        df_results = pd.DataFrame({
            'spearman_correlation': [spearman_corr],
            'pearson_correlation': [pearson_corr],
            'rmsd': [rmsd],
            'mantel_correlation': [mantel_corr],
            'mantel_p_value': [mantel_p_value]
        })
        return df_results
    
    
    def _convert_df_to_distance_matrix(self, df_with_scores: pd.DataFrame, score_col: str) -> DistanceMatrix:
        """
        Converts a DataFrame with pairwise scores into a distance matrix.
        Args:
            df (pd.DataFrame): DataFrame containing columns [PDBchain1, PDBchain2, score]
            score_col (str): The name of the column containing the scores.
        Returns:
            DistanceMatrix: A DistanceMatrix object representing the distance matrix.
        """
        structures = sorted(set(df_with_scores["PDBchain1"]).union(df_with_scores["PDBchain2"]))
        n = len(structures)
        df = pd.DataFrame(np.zeros((n, n)), index=structures, columns=structures)
        for _, row in df_with_scores.iterrows():
            i = row["PDBchain1"]
            j = row["PDBchain2"]
            df.loc[i, j] = row[score_col]
            df.loc[j, i] = row[score_col]
        np.fill_diagonal(df.values, 0)
        mat = DistanceMatrix(df.to_numpy())
        return mat
    
    
    def _compare_graphs(self, G_1: nx.Graph, G_2: nx.Graph) -> pd.DataFrame:
        """
        Compares 2 given graphs.
        Returns a dataframe with the comparison results.
        Args:
            G_1 (networkx.Graph): The first graph to compare.
            G_2 (networkx.Graph): The second graph to compare.
        Returns:
            list: A list containing dataframes with the comparison results for each graph.
        """
        results = []
        for G in [G_1, G_2]:
            metrics = self._compute_graph_metrics(G)
            results.append(metrics)
        return pd.DataFrame(results)
            
    
    def _compute_graph_metrics(self, G: nx.Graph) -> dict:
        """
        Computes various graph metrics for the given graph G.
        Returns a dictionary with the computed metrics.
        Args:
            G (networkx.Graph): The input graph.
        Returns:
            dict: A dictionary containing various graph metrics.
        """
        deg_per_node = dict(nx.degree(G))
        clustering_per_node = nx.clustering(G)
        betweenness_per_node = nx.betweenness_centrality(G)
        connected_components_lst = list(nx.connected_components(G))
        metrics = {
            "num_nodes": G.number_of_nodes(),
            "num_edges": G.number_of_edges(),
            "avg_degree": sum(deg_per_node.values()) / len(deg_per_node),
            "max_degree": max(deg_per_node.values()),
            "avg_clustering": sum(clustering_per_node.values()) / len(clustering_per_node),
            "avg_betweenness": sum(betweenness_per_node.values()) / len(betweenness_per_node),
            "num_connected_components": nx.number_connected_components(G),
            "average_shortest_path_length": nx.average_shortest_path_length(G)
        }
        return metrics
    

    def _build_graph_from_scores(self, df, coords_df):
        """
        Builds a 3D graph (NetworkX) from pairwise score data and node coordinates.

        Args:
            df (pd.DataFrame): DataFrame containing columns f.e. [PDBchain1, PDBchain2, score]
            coords_df (pd.DataFrame): DataFrame containing columns [PDBchain1, x, y, z]

        Returns:
            networkx.Graph: Graph where nodes represent PDB chains and edges represent pairwise scores.
                            Each node has a 'pos' attribute with 3D coordinates (x, y, z).
        """
        G = nx.Graph()
        for _, row in df.iterrows():
            G.add_edge(
                row.iloc[0], # node_x
                row.iloc[1], # node_y
                weight=row.iloc[2] # edge_length
            )
        coords_dict = dict(zip(coords_df.iloc[:, 0], 
                           zip(coords_df.iloc[:, 1], coords_df.iloc[:, 2], coords_df.iloc[:, 3])))
        nx.set_node_attributes(G, coords_dict, name="pos")
        return G
    
    
    def _get_cluster_labels(self, coordinates_df: pd.DataFrame) -> pd.DataFrame:
        """
        Finds cluster labels with the coordiantes_df containing the coordinates for the structure and sequence-based clans file.
        This method uses DBSCAN to find labels.
        Args:
            coordinates_df (pd.DataFrame): A DataFrame containing normalized x, y and z coordinates for each structure for the structure and sequence-based clans file.
        Returns:
            pd.DataFrame: A DataFrame containing the cluster labels for each structure for the structure and sequence-based clans file.
        """
        coords_struct = coordinates_df[["x_struct_z-scores", "y_struct_z-scores", "z_struct_z-scores"]]
        coords_seq = coordinates_df[["x_seq_z-scores", "y_seq_z-scores", "z_seq_z-scores"]]
        db_struct = DBSCAN(eps=10.0, min_samples=2).fit(coords_struct)
        db_seq = DBSCAN(eps=10.0, min_samples=2).fit(coords_seq)
        coordinates_df["cluster_label_struct"] = db_struct.labels_
        coordinates_df["cluster_label_seq"] = db_seq.labels_
        return coordinates_df
    
    
    def _map_cluster_label_to_structures(self, structures_with_labels: pd.DataFrame, label_col: str, structure_col: str) -> dict:
        """
        Maps cluster labels to the structures belonging to each cluster.
        Args:
            structures_with_labels (pd.DataFrame): A DataFrame containing the cluster labels for each structure.
                                                    columns = [structure, cluster_label]
            label_col (str): The name of the column containing the cluster labels.
            structure_col (str): The name of the column containing the structure identifiers.
        Returns:
            dict: A dictionary mapping each cluster label to a set of structures belonging to that cluster.
        """
        label_to_structures = {}
        labels = structures_with_labels[label_col].unique()
        labels = [l for l in labels if l != -1]  # exclude noise label (-1)
        for label in labels:
            structures = set(structures_with_labels[structures_with_labels[label_col] == label][structure_col])
            label_to_structures[label] = structures
        return label_to_structures
    
    
    def _compare_clusters(self, scores_combined: pd.DataFrame, structures_with_labels: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Compares clusters found with structure-based coordinates and sequence-based coordinates.

        Args:
            scores_combined (pd.DataFrame): A DataFrame containing pairwise scores with columns = ["PDBchain1", "PDBchain2", "score_struct", "score_seq"].
            structures_with_labels (pd.DataFrame): A DataFrame containing the cluster labels for each structure with columns = ["PDBchain1", "cluster_label_struct", "cluster_label_seq"].

        Returns:
            tuple[pd.DataFrame, pd.DataFrame]: A tuple containing two DataFrames: metrics per cluster and cluster overlap.
        """
        struct_scores = scores_combined[["PDBchain1", "PDBchain2", "score_struct"]].rename(columns={"score_struct": "score"})
        seq_scores = scores_combined[["PDBchain1", "PDBchain2", "score_seq"]].rename(columns={"score_seq": "score"})
        label_to_structures_struct = self._map_cluster_label_to_structures(structures_with_labels, "cluster_label_struct", "PDBchain1")
        label_to_structures_seq = self._map_cluster_label_to_structures(structures_with_labels, "cluster_label_seq", "PDBchain1")
        df_cluster_overlap = self._compute_cluster_overlap(label_to_structures_struct, label_to_structures_seq)
        df_cluster_metrics_struct = self._compute_metrics_per_cluster(struct_scores, label_to_structures_struct)
        df_cluster_metrics_seq = self._compute_metrics_per_cluster(seq_scores, label_to_structures_seq)
        df_cluster_metrics_struct["cluster_label"] = (df_cluster_metrics_struct["cluster_label"].astype(str) + "_struct")
        df_cluster_metrics_seq["cluster_label"] = (df_cluster_metrics_seq["cluster_label"].astype(str) + "_seq")
        df_metrics_per_cluster = pd.concat(
            [df_cluster_metrics_struct, df_cluster_metrics_seq],
            axis=0,
            ignore_index=True
        )
        return df_metrics_per_cluster, df_cluster_overlap
    
    
    def _compute_metrics_per_cluster(self, pairwise_scores: pd.DataFrame, label_to_structures: dict) -> pd.DataFrame:
        """
        Computes metrics for each cluster based on the pairwise_scores DataFrame.
        Args:
            pairwise_scores (pd.DataFrame): A DataFrame containing pairwise scores with columns = ["PDBchain1", "PDBchain2", "score"].
            label_to_structures (dict): A dictionary mapping cluster labels to sets of structures.
        Returns:
            pd.DataFrame: A DataFrame containing metrics for each cluster.
        """
        cluster_metrics = []
        for label, structures in label_to_structures.items():
            # get possible pairs within the cluster
            scores_of_cluster = pairwise_scores[pairwise_scores["PDBchain1"].isin(structures) & pairwise_scores["PDBchain2"].isin(structures)]
            cluster_metrics.append({
                "cluster_label": label,
                "num_structures": len(structures),
                "mean_score": scores_of_cluster["score"].mean(),
                "median_score": scores_of_cluster["score"].median(),
                "std_score": scores_of_cluster["score"].std()
            })
        return pd.DataFrame(cluster_metrics)
    
    
    def _compute_cluster_overlap(self, label_to_structures_struct: dict, label_to_structures_seq: dict) -> pd.DataFrame:
        """
        Computes the overlap between clusters found with structure-based coordinates and sequence-based coordinates.
        Args:
            structures_with_labels (pd.DataFrame): A DataFrame containing the cluster labels for each structure for the structural and sequence-based clans file.
                                                    columns = [structure, cluster_label_struct, cluster_label_seq]
        Returns:
            pd.DataFrame: A DataFrame containing the overlap metrics between structure-based and sequence-based clusters
            [cluster_label_struct, cluster_label_seq, num_structures_struct, num_structures_seq, num_overlap, num_union, jaccard_index].
        """
        overlap_results = []
        for label_struct, structures_struct in label_to_structures_struct.items():
            for label_seq, structures_seq in label_to_structures_seq.items():
                intersection = structures_struct & structures_seq
                union = structures_struct | structures_seq
                jacc = len(intersection) / len(union) if len(union) else 0
                overlap_results.append({
                    "cluster_label_struct": label_struct,
                    "cluster_label_seq": label_seq,
                    "num_structures_struct": len(structures_struct),
                    "num_structures_seq": len(structures_seq),
                    "num_overlap": len(intersection),
                    "num_union": len(union),
                    "jaccard_index": jacc
                })
        return pd.DataFrame(overlap_results) 
    
    
    def _get_combined_scores_and_coordinates_df(self, struct_clans_file: ClansFile, seq_clans_file: ClansFile):
        """
        Combines the scores of the structural clans file and the scores of the sequence clans file in one DataFrame.
        Args:
            struct_clans_file: A ClansFile object containing the structural scores.
            seq_clans_file: A ClansFile object containing the sequence scores.
        Returns:
            A list of dataframes. One with the merged scores with the first column for structural scores and the second column for sequence scores.
            And one for the cooridinates with each row representing the cooridantes of each PDBchain of the dataset.
        """
        structural_scores = struct_clans_file.scores
        structural_coordinates = struct_clans_file.get_coordinates()
        sequence_scores = seq_clans_file.scores
        sequence_coordinates = seq_clans_file.get_coordinates()
        merged_coordinates = pd.merge(structural_coordinates, sequence_coordinates, on=["PDBchain1"], suffixes=('_struct', '_seq'))
        merged_scores = pd.merge(structural_scores, sequence_scores, on=['PDBchain1', 'PDBchain2'], suffixes=('_struct', '_seq'))
        merged_scores = merged_scores[["PDBchain1", "PDBchain2", "score_struct", "score_seq"]].apply(pd.to_numeric, errors='coerce')
        return merged_scores, merged_coordinates
    
    
    def _match_clans_files_for_comparison(self, structural_clans_dir, sequence_clans_dir):
        """
        This method matches the clans files generated with structure similarity scores to the clans files generated with sequence similarity scores.
        It returns a dictionary which maps the structural clans file path to the sequence clans file path.
        """
        structural_to_sequence = {}
        for struct_file in os.listdir(structural_clans_dir):
            if not struct_file.endswith("_out.clans"):
                continue
            seq_file = struct_file
            structural_to_sequence[os.path.join(structural_clans_dir, struct_file)] = os.path.join(sequence_clans_dir, seq_file)
        return structural_to_sequence


    def _generate_input_output_files_dict(self, input_files_dir, output_files_dir):
        """
        Generates a dictionary mapping input files to output files.
        It appends "_out" to the output file names.
        Args:
            input_files_dir: Directory containing the input files.
            output_files_dir: Directory where the output files will be saved.
        Returns:
            dict: A dictionary mapping absolute input file paths to absolute output file paths.
        """
        input_output_files = {}

        for file in os.listdir(input_files_dir):
            # only use files that contain "cleaned" in the name
            if "cleaned" not in file:
                continue

            base_name = os.path.splitext(file)[0]
            output_file_name = f"{base_name}_out.clans"

            input_file_path = os.path.abspath(os.path.join(input_files_dir, file))
            output_file_path = os.path.abspath(os.path.join(output_files_dir, output_file_name))

            input_output_files[input_file_path] = output_file_path

        return input_output_files


    def _compute_struct_clans_files(self, scores_for_each_dataset: dict, fasta_files: dict):
        """
        Generates clans files for each fasta file based on the computed scores.

        Args:
            scores_for_each_dataset (dict): A dictionary containing the scores for each dataset.
            fasta_files: A dictionary mapping dataset names to paths of the cleaned fasta files.
        """
        for dataset_name, cleaned_dataset_path in fasta_files.items():
            self.struct_clans_generator.generate_clans_file(scores_for_each_dataset[dataset_name], cleaned_dataset_path)


    def _compute_struct_scores(self, structures_dir: str) -> dict:
        """
        Computes pairwise structure similarity scores for the given structures.
        
        Args:
            structures_dir: A directory with directories of pdbs for each dataset
        Returns:
            dict: A dictionary containing the scores for each dataset.
        """
        scores = {}
        for sub_dir in os.listdir(structures_dir):
            path_to_subdir = os.path.join(structures_dir, sub_dir)
            print(f"Computing scores for dataset {sub_dir}...")
            scores[sub_dir] = self.scores_computer.run(self.tool_type, path_to_subdir)
        return scores
    
    
    def _generate_fasta_from_uids_with_regions_for_each_dataset(self, uids_with_regions_for_each_dataset: dict, datasets_dir: str, datasets_file_type: InputFileType) -> dict:
        """
        Generates fasta files from the given uids with regions for each dataset.
        Args:
            uids_with_regions_for_each_dataset: A dictionary containing uids with regions for each dataset.
            datasets_dir: Directory where the datasets are stored.
            datasets_file_type: Specifies the type of the dataset files.
        Returns:
            dict: A dictionary mapping dataset names to paths of the generated fasta files.
        """
        paths_to_cleaned_datasets = {}
        for dataset in os.listdir(datasets_dir):
            dataset_name = os.path.splitext(dataset)[0]
            cleaned_dataset_path = os.path.join(datasets_dir, f"{dataset_name}_cleaned.fasta")
            uids_with_regions = uids_with_regions_for_each_dataset[dataset_name]
            if datasets_file_type == InputFileType.FASTA:
                original_fasta_path = os.path.join(datasets_dir, dataset)
            else:
                original_fasta_path = None
            cleaned_dataset_path = generate_fasta_from_uids_with_regions(uids_with_regions, cleaned_dataset_path, original_fasta_path)
            paths_to_cleaned_datasets[dataset_name] = cleaned_dataset_path
        return paths_to_cleaned_datasets


    def _download_structures(self, path_to_datasets: str, out_dir: str, datasets_file_type: InputFileType) -> dict:
        """
        Downloads structure files for the datasets.
        For each dataset, the corresponding structure files are downloaded in a separated directory in the self.structures_dir.
        Args:
            path_to_datasets: path to the datasets for which to download the corresponding structure files
            out_dir: path to the directory where the structure files are saved
        Returns:
            dict: Containing downloaded uids together with their regions for each dataset.
        """
        uids_with_regions_for_each_dataset = {}
        print("Downloading structure files for the datasets...")
        datasets = [f for f in os.listdir(path_to_datasets)]
        for dataset in datasets:
            dataset_path = os.path.join(path_to_datasets, dataset)
            # create dir for each dataset
            dataset_name = os.path.splitext(dataset)[0]
            dir_for_dataset_structures = os.path.join(out_dir, dataset_name)
            reset_dir_content(dir_for_dataset_structures)
            print(f"Downloading structure files for dataset {dataset}...")
            uids_with_regions = fetch_pdbs(dataset_path, datasets_file_type, dir_for_dataset_structures)
            uids_with_regions_for_each_dataset[dataset_name] = uids_with_regions
        return uids_with_regions_for_each_dataset
    