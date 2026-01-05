import sys
import os
sys.path.append(os.path.abspath(".."))
from ToolType import ToolType
from InputFileType import InputFileType
from main import create_clans_file
from recovered_CLANS.utils_old_clans import generate_clans_file_seq_based
from ClansFileGenerator import ClansFileGenerator


class ScoresEvaluator:
    def __init__(self, working_dir: str, data: str, input_file_type: InputFileType):
        self.working_dir = working_dir
        self.blast_dir = os.path.join(working_dir, "blast_temp")
        self.data = data
        self.input_file_type = input_file_type
        self.clans_file_generator = ClansFileGenerator(self.working_dir)


    def generate_clans_files(self, tool: ToolType, score: str | None) -> tuple[str, str]:
        struct_clans_file_path, cleaned_input_file_path = create_clans_file(self.data, self.input_file_type, tool, score)
        seq_clans_file_path = generate_clans_file_seq_based(cleaned_input_file_path, self.blast_dir, self.clans_file_generator, self.working_dir)
        return struct_clans_file_path, seq_clans_file_path