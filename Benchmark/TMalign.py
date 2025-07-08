from Tool import Tool
import os
import subprocess
"""
This class extends the Tool class to implement the TMalign tool for protein structure alignment.
"""


class TMalign(Tool):
    def __init__(self, name, description, command):
        super().__init__(name, description, command)
