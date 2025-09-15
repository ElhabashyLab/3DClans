from Dataset_Generator.DatasetGenerator import DatasetGenerator
from fasta2PDB import *
from StructSimComputer import StructSimComputer
from ClansFileGenerator import ClansFileGenerator
from ToolType import ToolType
import os


class ScoresEvaluator:
    """
    This class is responsible for evaluating the structure similarity scores by comparing them to sequence similarity scores.
    The datasets used for the evaluation are generated with the DatasetGenerator class with the given seeds.
    With the generated data, the structure similarity scores and sequence similarity scores are computed and clans files anre generated.
    The comparison is done with the scores as well as with inferred graphs and clusters.
    """
    def __init__(self, number_of_datasets, size_of_datasets, seeds):
        self.number_of_datasets = number_of_datasets
        self.size_of_datasets = size_of_datasets
        self.seeds = seeds
        self.working_dir = "Scores_Evaluation"
        self._set_up_dirs()
        self.dataset_generator = DatasetGenerator(self.generated_datasets_dir)
        self.scores_computer = StructSimComputer()
        self.tool_type = ToolType.FOLDSEEK
        self.clans_generator = ClansFileGenerator(self.clans_files_dir)      


    def _set_up_dirs(self):
        """
        Sets up the necessary directories for the evaluation.
        Creates directories for generated datasets, PDB files and clans files and overwrites them if they already exist.
        """
        self.generated_datasets_dir = f"{self.working_dir}/generated_datasets"
        self.pdbs_dir = f"{self.working_dir}/PDBs"
        self.clans_files_dir = f"{self.working_dir}/clans_files"
        for dir_path in [self.generated_datasets_dir, self.pdbs_dir, self.clans_files_dir]:
            delete_dir_content(dir_path)
        os.makedirs(self.pdbs_dir, exist_ok=True)
        os.makedirs(self.clans_files_dir, exist_ok=True)
        os.makedirs(self.generated_datasets_dir, exist_ok=True)


    def run_evaluation(self):
        """
        Runs the evaluation by generating datasets, computing scores, generating clans files and comparing the results.
        """
        self._check_for_correct_input()
        self._generate_datasets()
        self._download_pdbs()
        scores_for_each_dataset = self._compute_scores()
        self._compute_clans_files(scores_for_each_dataset)
        
        
    def _compute_clans_files(self, scores_for_each_dataset: dict):
        """
        Generates clans files for each dataset based on the computed scores.

        Args:
            scores_for_each_dataset (dict): A dictionary containing the scores for each dataset.
        """
        for i in range(self.number_of_datasets):
            cleaned_dataset_i_path = os.path.join(self.generated_datasets_dir, f"dataset_{i+1}_cleaned.fasta")
            self.clans_generator.generate_clans_file(scores_for_each_dataset[i], cleaned_dataset_i_path)

    
    def _compute_scores(self):
        """
        Computes the structure similarity scores for the generated datasets.
        
        Returns:
            dict: A dictionary containing the scores for each dataset.
        """
        scores = {}
        for i in range(self.number_of_datasets):
            dataset_path = os.path.join(self.generated_datasets_dir, f"dataset_{i+1}.fasta")
            pdb_dir_for_dataset = os.path.join(self.pdbs_dir, f"dataset_{i+1}")
            print(f"Computing scores for dataset {i+1}/{self.number_of_datasets}...")
            scores[i] = self.scores_computer.run(self.tool_type, pdb_dir_for_dataset)
        return scores


    def _download_pdbs(self):
        """
        Downloads the PDB files for the generated datasets.
        For each dataset, the corresponding PDB files are downloaded in a separated directory in the self.pdbs_dir.
        """
        print("Downloading PDB files for the generated datasets...")
        dataset_files = [f for f in os.listdir(self.generated_datasets_dir) if f.endswith(".fasta")]
        for dataset in dataset_files:
            dataset_path = os.path.join(self.generated_datasets_dir, dataset)
            # create pdb dir for each dataset
            pdb_dir_for_dataset = os.path.join(self.pdbs_dir, dataset.replace(".fasta", ""))
            if not os.path.exists(pdb_dir_for_dataset):
                os.makedirs(pdb_dir_for_dataset)
            print(f"Downloading PDB files for dataset {dataset}...")
            fetch_pdbs(dataset_path, pdb_dir_for_dataset)
    
    
    def _generate_datasets(self):
        """
        Generates the datasets (fasta_files) by calling the dataset generator.
        The datasets are generated in the self.generated_datasets_dir.
        """
        n_cl_per_dataset = len(self.seeds) // self.number_of_datasets
        print(f"Generating {self.number_of_datasets} datasets with {self.size_of_datasets} sequences each...")
        for i in range(self.number_of_datasets):
            print(f"Generating dataset {i+1}/{self.number_of_datasets}...")
            seeds_for_current_dataset = self.seeds[i*n_cl_per_dataset:(i+1)*n_cl_per_dataset]
            self.dataset_generator.generate(self.size_of_datasets, n_cl_per_dataset, seeds_for_current_dataset, f"dataset_{i+1}.fasta")


    def _check_for_correct_input(self):
        """Checks if the input parameters are correct.

        Raises:
            ValueError: If number_of_datasets is less than or equal to 0.
            ValueError: If the number of seeds is less than the number of datasets.

        Returns:
            None
        """
        if self.number_of_datasets <= 0:
            raise ValueError("Number of datasets must be greater than 0.")
        if len(self.seeds) < self.number_of_datasets:
            raise ValueError("Number of seeds must be greater than or equal to number of datasets.")
    
    
# test
evaluator = ScoresEvaluator(2, 20, ["Q99895", "P42212"])
evaluator.run_evaluation()
