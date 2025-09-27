from Dataset_Generator.DatasetGenerator import DatasetGenerator
from old_clans.utils_old_clans import run_multiple_clans_headless
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
    def __init__(self, number_of_datasets, size_of_datasets, seeds, path_to_recovered_clans: str = "/home/aronw/Development/clans-recovered", use_existing_datasets: bool = False):
        self._check_for_correct_input(number_of_datasets, seeds)
        self.size_of_datasets = size_of_datasets
        self.working_dir = "Scores_Evaluation"
        self._set_up_dirs()
        self.dataset_generator = DatasetGenerator(self.generated_datasets_dir)
        self.scores_computer = StructSimComputer()
        self.tool_type = ToolType.FOLDSEEK
        self.clans_generator = ClansFileGenerator(self.clans_files_dir)   
        self.path_to_recovered_clans = path_to_recovered_clans
        
        
    def _set_up_dirs(self):
        """
        Sets up the necessary directories for the evaluation.
        Creates directories for generated datasets, PDB files and clans files and overwrites them if they already exist.
        """
        self.generated_datasets_dir = f"{self.working_dir}/generated_datasets"
        self.pdbs_dir = f"{self.working_dir}/PDBs"
        self.clans_files_dir = f"{self.working_dir}/clans_files_structsim" # dir for output of clans files generated with structure similarity scores
        self.clans_files_seqsim_dir = f"{self.working_dir}/clans_files_seqsim" # dir for output of clans runs with old software
        for dir_path in [self.generated_datasets_dir, self.pdbs_dir, self.clans_files_dir]:
            delete_dir_content(dir_path)


    def _initialize_evaluation(self):
        """
        Generates datasets (fasta_files), downloads the corresponding PDB files, computes structure similarity scores and creates clans files from those.
        After it runs the old recovered clans on the generated clans files and saves the updated clans files.
        This function should be called before running the actual evaluation.
        """
        print("Running scores evaluation with generated datasets ...")
        self._generate_datasets(self.number_of_datasets, self.size_of_datasets, self.seeds)
        self._download_pdbs(self.generated_datasets_dir, self.pdbs_dir)
        scores_for_each_dataset = self._compute_scores(self.pdbs_dir, self.number_of_datasets)
        self._compute_clans_files(scores_for_each_dataset, self.generated_datasets_dir, self.number_of_datasets)
        # input clans files will be overwritten with updated clans files
        input_output_files_for_struct_clans_files = self.generate_input_output_files_dict(self.clans_files_dir)
        run_multiple_clans_headless(self.path_to_recovered_clans, input_output_files_for_struct_clans_files)
        
    
    def generate_input_output_files_dict(self, input_files_dir, overwrite: bool = False):
        """
        Generates a dictionary mapping input files to output files.
        If overwrite is True, the output files will be saved in the same directory as the input files and overwrite them.
        If overwrite is False, the output files will be saved with a "_out" suffix.
        Args:
            input_files_dir: Directory containing the input files.
            overwrite: Whether to overwrite the input files or save the output files with a "_out" suffix.
        Returns:
            dict: A dictionary mapping input files to output files.
        """
        input_output_files = {}
        for file in os.listdir(input_files_dir):
            if overwrite:
                input_output_files[os.path.abspath(os.path.join(input_files_dir, file))] = os.path.abspath(os.path.join(input_files_dir, file))
            else:
                input_output_files[os.path.abspath(os.path.join(input_files_dir, file))] = os.path.abspath(os.path.join(input_files_dir, file.replace(".clans", "_out.clans")))
        print(f"Generated input-output files dict: {input_output_files}")
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
            self.clans_generator.generate_clans_file(scores_for_each_dataset[i], cleaned_dataset_i_path)

    
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
            
    
# test
evaluator = ScoresEvaluator(2, 20, ["Q99895", "P42212"])
evaluator._initialize_evaluation()
