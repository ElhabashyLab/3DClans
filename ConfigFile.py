import shlex
from InputFileType import InputFileType


class ConfigFile():
    def __init__(self, filepath: str):
        self.filepath = filepath
        self.InputFileType = InputFileType.CONF
        self.documentation = "# Configuration file for CLANS\n# Each line contains a key-value pair in the format -<key> <value>.\n# Lines starting with '#' are comments and will be ignored."


    def write_config(self, arguments: dict):
        """
        Writes the configuration file with the given arguments to the specified filepath.

        Args:
            arguments (dict): A dictionary containing the configuration key-value pairs.
        Returns: None
        """
        with open(self.filepath, 'w') as f:
            f.write(self.documentation + "\n")
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
    
    
# tests
#test_conf = ConfigFile("test.conf")
#test_conf.write_config({"Arg1": "Value1", "Arg2": "Value2"})
#config_data = test_conf.read_config()
#print(config_data)  # Expected output: {'Arg1': 'Value1', 'Arg2': 'Value2'}