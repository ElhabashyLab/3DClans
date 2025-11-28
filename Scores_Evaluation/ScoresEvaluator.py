import sys
import os
sys.path.append(os.path.abspath(".."))
from Dataset_Generator.DatasetGenerator import DatasetGenerator
from old_clans.utils_old_clans import run_clans_headless
from utils_for_PDB import copy_dir_content, reset_dir_content
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
        self.scores_computer = StructSimComputer()
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
        uids_with_regions_for_each_dataset = self._download_pdbs(self.datasets_dir, self.structures_dir, datasets_file_type)
        paths_to_cleaned_datasets = self._generate_fasta_from_uids_with_regions_for_each_dataset(uids_with_regions_for_each_dataset, self.datasets_dir, datasets_file_type)
        scores_for_each_dataset = self._compute_struct_scores(self.structures_dir)
        self._compute_struct_clans_files(scores_for_each_dataset, paths_to_cleaned_datasets)
        # clustering clans files with recovered clans
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
        Evaluates the clustered clans files generated with structure similarity scores and sequence similarity scores.
        Returns a dataframe with the evaluation results consisting of numerical comparison, graph comparison and cluster comparison.
        Args:
            structural_to_sequence_clans_files: A dictionary mapping structural clans file paths to sequence clans file paths.
        """
        df_evaluation_results = []
        df_numerical_comparisons = pd.DataFrame()
        figures = []
        for struct_clans_file, seq_clans_file in structural_to_sequence_clans_files.items():
            print(f"Evaluating structural clans file {struct_clans_file} with sequence clans file {seq_clans_file}...")
            parsed_struct_clans_file = ClansFileGenerator.parse_clans_file(struct_clans_file)
            parsed_seq_clans_file = ClansFileGenerator.parse_clans_file(seq_clans_file)
            merged_scores_df, merged_coordinates_df = self._prepare_scores(parsed_struct_clans_file, parsed_seq_clans_file)
            print(merged_scores_df, merged_coordinates_df)
            
            # numerical evaluation
            figures.append(self._plot_similarity_scatter(merged_scores_df, x_col="score_struct_log10_min-max-scaled", y_col="score_seq_log10_min-max-scaled", prevent_display=False))
            df_numerical_comparison = self._compare_numerically(merged_scores_df, ["score_struct_log10", "score_seq_log10"])
            df_numerical_comparisons = pd.concat([df_numerical_comparisons, df_numerical_comparison], ignore_index=True)
            
            # graph evaluation
            coordinates_for_struct = merged_coordinates_df[["PDBchain1", "x_struct", "y_struct", "z_struct"]]
            coordinates_for_seq = merged_coordinates_df[["PDBchain1", "x_seq", "y_seq", "z_seq"]]
            G_struct = self._build_graph_from_scores(merged_scores_df[["PDBchain1", "PDBchain2", "score_struct"]], coordinates_for_struct)
            G_seq = self._build_graph_from_scores(merged_scores_df[["PDBchain1", "PDBchain2", "score_seq"]], coordinates_for_seq)
            #df_graph_comparison = self._compare_graphs(G_struct, G_seq)
            
        print("Evaluation completed.")
        return df_evaluation_results, figures
        
        
    def display_evaluation_results(self, df_numerical_comparisons, figures):
        print("Scatterplots: ")
        for fig in figures:
            fig
        print("Numerical comparison results:")
        print(df_numerical_comparisons)
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
            
    
    def _prepare_scores(self, struct_clans_file: ClansFile, seq_clans_file: ClansFile) -> tuple[pd.DataFrame, pd.DataFrame]:
        """
        Prepares the scores of parsed clans file for further analyis. On of the Clans files should contain structural scores and the other sequence-based scores.
        It merges the scores into one DataFrame,
        log-transforms and normalizes the scores (min-max-scaling and z-scores are computed).

        Args:
            struct_clans_file (ClansFile): A ClansFile Object based on structural similarity.
            seq_clans_file (ClansFile): A ClansFIle Object based on sequence similarity.

        Returns:
            pd.DataFrame: A DataFrame containing the original scores of both ClansFiles,
            log-transformed scores, min-max scaled-scores and z-scores and a DataFrame containing the coordinates for each PDBchain.
        """
        merged_scores_df, merged_coordinates = self._get_merged_scores_and_cooridinates_df(struct_clans_file, seq_clans_file)
        merged_scores_df = self._log_transform(merged_scores_df, ["score_struct", "score_seq"])
        merged_scores_df = self._normalize_scores(merged_scores_df, ["score_struct_log10", "score_seq_log10"], "z-score")
        merged_scores_df = self._normalize_scores(merged_scores_df, ["score_struct_log10", "score_seq_log10"], "min-max")
        return merged_scores_df, merged_coordinates
        
        
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
        

    def _compare_numerically(self, df_combined: pd.DataFrame, columns: list) -> pd.DataFrame:
        """
        Compares the scores of 2 given columns numerically.
        Returns a dataframe with the comparison results.
        The scores of the given columns are compared with the metrics: 
        Correlation(Spearman and Pearson), RMSD
        Args:
            df_combined: A dataframe with the merged scores of the structural clans file and the sequence clans file.
            columns: A list of 2 columns contained in df_combined.
        Returns:
            df_results: A dataframe containg Correlation(Spearman and Pearson) and RMSD values for the scores of the given columns.
        """
        df_combined = df_combined[columns].apply(pd.to_numeric)
        spearman_corr = df_combined[columns[0]].corr(df_combined[columns[1]], method='spearman')
        pearson_corr = df_combined[columns[0]].corr(df_combined[columns[1]], method='pearson')
        rmsd = ((df_combined[columns[0]] - df_combined[columns[1]]) ** 2).mean() ** 0.5
        df_results = pd.DataFrame({
            'spearman_correlation': [spearman_corr],
            'pearson_correlation': [pearson_corr],
            'rmsd': [rmsd]
        })
        return df_results
    
    
    def _compare_graphs(self, G_1, G_2):
        """
        Compares the graphs inferred from the structural clans file with the graphs inferred from the sequence clans file.
        Returns a dataframe with the comparison results.
        """
        
        raise NotImplementedError("Graph comparison not implemented yet.")
    

    def _build_graph_from_scores(self, df, coords_df):
        """
        Builds a 3D graph (NetworkX) from pairwise score data and node coordinates.

        Args:
            df (pd.DataFrame): DataFrame containing columns f.e. [PDBchain1, PDBchain2, score]
            coords_df (pd.DataFrame): DataFrame containing columns [PDBchain, x, y, z]

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
    
    
    def _compare_clusters(self, struct_clans_file: ClansFile, seq_clans_file: ClansFile):
        """
        Compares the clusters inferred from the structural clans file with the clusters inferred from the sequence clans file.
        Returns a dataframe with the comparison results.
        """
        raise NotImplementedError("Cluster comparison not implemented yet.")
    
    
    def _get_merged_scores_and_cooridinates_df(self, struct_clans_file: ClansFile, seq_clans_file: ClansFile):
        """
        Merges the scores of the structural clans file with the scores of the sequence clans file.
        Returns a dataframe with the merged scores.
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


    def _download_pdbs(self, path_to_datasets: str, out_dir: str, datasets_file_type: InputFileType) -> dict:
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
    