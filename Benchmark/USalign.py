from Tool import Tool
"""
This class extends the Tool class to implement the USalign tool for protein structure alignment.
"""


class USalign(Tool):
    def __init__(self, name, description, command):
        super().__init__(name, description, command)
