import shlex
from clans3d.core.input_file_type import InputFileType


class ConfigFile():
    def __init__(self, filepath: str, documentation: str|None = None):
        self.filepath = filepath
        self.InputFileType = InputFileType.CONF
        self.documentation = documentation


    def write_config(self, arguments: dict):
        """
        Writes the configuration file with the given arguments to self.filepath.

        Args:
            arguments (dict): A dictionary containing the configuration key-value pairs.
        Returns: None
        """
        with open(self.filepath, 'w') as f:
            if self.documentation:
                f.write("# " + self.documentation + "\n")
            for key, value in arguments.items():
                f.write(f"-{key} {value}\n")


    def read_config(self) -> dict[str, str]:
        """
        Reads the configuration file and returns a dictionary of key-value pairs.
        Lines starting with '#' or empty lines are ignored.
        Inline comments are supported.
        Arguments are expected in the format -<key> <value>.

        Returns:
            dict: A dictionary containing the configuration key-value pairs.
        """
        config = {}
        with open(self.filepath, "r") as f:
            for lineno, line in enumerate(f, start=1):
                line = line.strip()
                if not line or line.startswith("#"): # ignore comments and empty lines
                    continue
                # get key
                tokens = shlex.split(line, comments=True)
                if not tokens:
                    continue
                key_token = tokens[0]
                if not key_token.startswith("-"):
                    raise ValueError(
                        f"Invalid config entry on line {lineno}: "
                        f"expected '-<key> <value>'"
                    )
                key = key_token.lstrip("-")
                if not key:
                    raise ValueError(
                        f"Invalid config entry on line {lineno}: empty key"
                    )
                # get value
                value = tokens[1] if len(tokens) > 1 else ""

                config[key] = value
        return config
