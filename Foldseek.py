from StructSimTool import StructSimTool
"""
This class extends the StructSimTool class to implement the Foldseek tool for protein structure comparison.
"""


class USalign(StructSimTool):
    def __init__(self):
        description = "A tool for protein structure comparison using Foldseek."
        super().__init__("foldseek", description)
        