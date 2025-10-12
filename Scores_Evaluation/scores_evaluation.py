from Dataset_Generator.DatasetGenerator import DatasetGenerator
from old_clans.utils_old_clans import run_clans_headless
from fasta2PDB import *
from StructSimComputer import StructSimComputer
from ClansFileGenerator import ClansFileGenerator
from ToolType import ToolType
import os
import pandas as pd 


CLUSTER_ROUNDS = 100000

class ScoresEvaluator:
    """
    This class is responsible for evaluating the structure similarity scores by comparing them to sequence similarity scores.
    The datasets used for the evaluation are generated with the DatasetGenerator class with the given seeds.
    With the generated data, the structure similarity scores and sequence similarity scores are computed and clans files anre generated.
    The comparison is done with the scores as well as with inferred graphs and clusters.
    """
    def __init__(self, number_of_datasets, size_of_datasets, seeds, path_to_recovered_clans: str = "/home/aronw/Development/clans-recovered"):
        self._check_for_correct_input(number_of_datasets, seeds)
        self.size_of_datasets = size_of_datasets
        self.working_dir = "Scores_Evaluation"
        self._set_up_dirs()
        self.dataset_generator = DatasetGenerator(self.generated_datasets_dir)
        self.scores_computer = StructSimComputer()
        self.tool_type = ToolType.FOLDSEEK
        self.struct_clans_generator = ClansFileGenerator(self.clans_files_structsim_dir)
        self.seq_clans_generator = ClansFileGenerator(self.clans_files_seqsim_dir)
        self.path_to_recovered_clans = path_to_recovered_clans
        
        
    def _set_up_dirs(self):
        """
        Sets up the necessary directories for the evaluation.
        Creates directories for generated datasets, pdbs, clans files and blast files.
        Deletes the content of the directories if they already exist.
        """
        abs_path_working_dir = os.path.abspath(self.working_dir)
        self.generated_datasets_dir = f"{abs_path_working_dir}/generated_datasets"
        self.pdbs_dir = f"{abs_path_working_dir}/PDBs"
        self.clans_files_structsim_dir = f"{abs_path_working_dir}/clans_files_structsim" # dir for clans files generated with structure similarity scores
        self.clans_files_seqsim_dir = f"{abs_path_working_dir}/clans_files_seqsim" # dir for clans files generated with sequence similarity scores
        self.blast_dir = f"{abs_path_working_dir}/blast_files"
        for dir_path in [self.generated_datasets_dir, self.pdbs_dir, self.clans_files_structsim_dir, self.blast_dir, self.clans_files_seqsim_dir]:
            delete_dir_content(dir_path)


    def _initialize_evaluation(self):
        """
        Generates datasets (fasta_files), downloads the corresponding PDB files, computes structure similarity scores and creates clans files from those.
        After it runs the old recovered clans on the generated clans files and saves the updated clans files.
        This function should be called before running the actual evaluation.
        """
        print("Initializing evaluation...")
        self._generate_datasets(self.number_of_datasets, self.size_of_datasets, self.seeds)
        self._download_pdbs(self.generated_datasets_dir, self.pdbs_dir)
        scores_for_each_dataset = self._compute_scores(self.pdbs_dir, self.number_of_datasets)
        self._compute_clans_files(scores_for_each_dataset, self.generated_datasets_dir, self.number_of_datasets)
        print("Running recovered clans.jar on the generated structural clans files...")
        input_output_dict_structural = self._generate_input_output_files_dict(self.clans_files_structsim_dir, self.clans_files_structsim_dir)
        run_clans_headless(self.path_to_recovered_clans, input_output_dict_structural, input_file_type="clans", rounds=CLUSTER_ROUNDS)
        self._remove_colorcutoffs_colorarr(self.clans_files_structsim_dir)
        print("Running recovered clans.jar with sequences similarity scores...")
        input_output_dict_sequence = self._generate_input_output_files_dict(self.generated_datasets_dir, self.clans_files_seqsim_dir)
        run_clans_headless(self.path_to_recovered_clans, input_output_dict_sequence, input_file_type="fasta", blast_dir=self.blast_dir, clans_generator=self.seq_clans_generator, rounds=CLUSTER_ROUNDS)
        self._remove_colorcutoffs_colorarr(self.clans_files_seqsim_dir)
        print("Evaluation initialized. Clans files generated and clustered for each dataset with structure similarity scores and sequence similarity scores.")


    def _remove_colorcutoffs_colorarr(self, clans_files_dir):
        """
        Fixes a temporary bug in the recovered clans code.
        This function is a temporary workaround and should be removed once the bug is fixed in the recovered clans code.
        After the clustering process the clans files contain 2 lines starting with colorcutoffs and colorarr with missing values.
        This function removes those lines from the clans files in the given directory so they can be loaded again.
        """
        for file in os.listdir(clans_files_dir):
            if file.endswith("_out.clans"):
                file_path = os.path.join(clans_files_dir, file)
                with open(file_path, "r") as f:
                    lines = f.readlines()
                with open(file_path, "w") as f:
                    for line in lines:
                        if line.startswith("colorcutoffs") or line.startswith("colorarr"):
                            continue # remove buggy lines
                        f.write(line)
    
    
    def evaluate(self, structural_to_sequence_clans_files: dict):
        """
        Evaluates the clustered clans files generated with structure similarity scores and sequence similarity scores.
        Returns a dataframe with the evaluation results consisting of numerical comparison, graph comparison and cluster comparison.
        Args:
            structural_to_sequence_clans_files: A dictionary mapping structural clans file paths to sequence clans file paths.
        """
        df_evaluation_combined = pd.DataFrame()
        for struct_clans_file, seq_clans_file in structural_to_sequence_clans_files.items():
            print(f"Evaluating structural clans file {struct_clans_file} with sequence clans file {seq_clans_file}...")
            df_numerical_comparison = self._compare_nmerically(struct_clans_file, seq_clans_file)
            df_graphs_comparison = self._compare_graphs(struct_clans_file, seq_clans_file)
            df_clusters_comparison = self._compare_clusters(struct_clans_file, seq_clans_file)
            df_evaluation_combined = pd.concat([df_numerical_comparison, df_graphs_comparison, df_clusters_comparison], axis=1)
        print("Evaluation completed.")
        return df_evaluation_combined
    
    
    def _compare_nmerically(self, struct_clans_file, seq_clans_file):
        """
        Compares the structural clans file with the sequence clans file numerically.
        Returns a dataframe with the comparison results.
        """
        raise NotImplementedError("Numerical comparison not implemented yet.")
    
    
    def _compare_graphs(self, struct_clans_file, seq_clans_file):
        """
        Compares the graphs inferred from the structural clans file with the graphs inferred from the sequence clans file.
        Returns a dataframe with the comparison results.
        """
        raise NotImplementedError("Graph comparison not implemented yet.")
    
    
    def _compare_clusters(self, struct_clans_file, seq_clans_file):
        """
        Compares the clusters inferred from the structural clans file with the clusters inferred from the sequence clans file.
        Returns a dataframe with the comparison results.
        """
        raise NotImplementedError("Cluster comparison not implemented yet.")
    
    
    def match_clans_files_for_comparison(self, structural_clans_dir, sequence_clans_dir):
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


    def _compute_clans_files(self, scores_for_each_dataset: dict, generated_datasets_dir, number_of_datasets):
        """
        Generates clans files for each dataset based on the computed scores.

        Args:
            scores_for_each_dataset (dict): A dictionary containing the scores for each dataset.
            generated_datasets_dir: A path to the generated datasets
            number_of_datasets: number of generated datasets
        """
        for i in range(number_of_datasets):
            cleaned_dataset_i_path = os.path.join(generated_datasets_dir, f"dataset_{i+1}_cleaned.fasta")
            self.struct_clans_generator.generate_clans_file(scores_for_each_dataset[i], cleaned_dataset_i_path)

    
    def _compute_scores(self, pdbs, number_of_datasets):
        """
        Computes the structure similarity scores for the generated datasets.
        
        Args:
            pdbs: A directory with directorys of pdbs for each dataset
        Returns:
            dict: A dictionary containing the scores for each dataset.
        """
        scores = {}
        for i in range(number_of_datasets):
            pdb_dir_for_dataset = os.path.join(pdbs, f"dataset_{i+1}")
            print(f"Computing scores for dataset {i+1}/{number_of_datasets}...")
            scores[i] = self.scores_computer.run(self.tool_type, pdb_dir_for_dataset)
        return scores


    def _download_pdbs(self, generated_datasets_dir, out_dir):
        """
        Downloads the PDB files for the generated datasets.
        For each dataset, the corresponding PDB files are downloaded in a separated directory in the self.pdbs_dir.
        Args:
            path_to_datasets: path to the datasets for which to download the corresponding pdbs
            out_dir: path to the directory where the pdbs are saved
        """
        print("Downloading PDB files for the generated datasets...")
        dataset_files = [f for f in os.listdir(generated_datasets_dir) if f.endswith(".fasta")]
        for i in range(len(dataset_files)):
            dataset_path = os.path.join(generated_datasets_dir,f"dataset_{i+1}.fasta")
            # create pdb dir for each dataset
            pdb_dir_for_dataset = os.path.join(out_dir, f"dataset_{i+1}")
            delete_dir_content(pdb_dir_for_dataset)
            print(f"Downloading PDB files for dataset dataset_{i+1}...")
            fetch_pdbs(dataset_path, pdb_dir_for_dataset)
    
    
    def _generate_datasets(self, n, size, seeds):
        """
        Generates n datasets (fasta_files) by calling the dataset generator.
        The datasets are saved in the out_dir.
        
        Args
            n: number of datasets to generate
            size: size of each dataset
            seeds: seed sequences for the Dataset-Generator
            out_dir: directory where to save the datasets
        """
        n_cl_per_dataset = len(seeds) // n
        print(f"Generating {n} datasets with {size} sequences each...")
        for i in range(n):
            print(f"Generating dataset {i+1}/{n}...")
            seeds_for_current_dataset = seeds[i*n_cl_per_dataset:(i+1)*n_cl_per_dataset]
            self.dataset_generator.generate(size, n_cl_per_dataset, seeds_for_current_dataset, f"dataset_{i+1}.fasta")


    def _check_for_correct_input(self, number_of_datasets, seeds):
        """Checks if the input parameters are correct.

        Raises:
            ValueError: If number_of_datasets is less than or equal to 0.
            ValueError: If the number of seeds is less than the number of datasets.

        Returns:
            None
        """
        if number_of_datasets <= 0:
            raise ValueError("Number of datasets must be greater than 0.")
        elif len(seeds) < number_of_datasets:
            raise ValueError("Number of seeds must be greater than or equal to number of datasets.")
        else:
            self.number_of_datasets = number_of_datasets
            self.seeds = seeds
            
    
# test ["P68871", "Q99895", "P42212", "P00734", "P69905", "P0A6F5"]
# evaluator = ScoresEvaluator(3, 30,  ["P68871", "Q99895", "P42212", "P00734", "P69905", "P0A6F5"])
# evaluator._initialize_evaluation()
# structural_to_sequence = evaluator.match_clans_files_for_comparison("Scores_Evaluation/clans_files_structsim", "Scores_Evaluation/clans_files_seqsim")
