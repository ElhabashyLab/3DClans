"""
This module defines a base class for tools that can be used in a benchmarking context.
"""


class Tool:
    def __init__(self, name: str, description: str, command: str):
        self.name = name
        self.description = description
        self.command = command

    def run(self):
        """
        Runs the tool using the specified command.
        This method should be overridden by subclasses to implement specific tool functionality.
        """
        raise NotImplementedError("Subclasses should implement this method.")



