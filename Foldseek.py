from Benchmark.StructSimTool import StructSimTool
"""
This class extends the Tool class to implement the Foldseek tool for protein structure comparison.
"""


class TMalign(StructSimTool):
    def __init__(self, name, description, command):
        super().__init__(name, description, command)