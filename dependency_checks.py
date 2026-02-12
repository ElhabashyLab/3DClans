import shutil
from typing import Optional
from ToolType import ToolType


def check_external_tool(tool_name: str) -> Optional[str]:
    """
    Check if an external tool is available in PATH.
    
    Args: 
        tool_name (str): The name of the tool executable to check.
        
    Returns:
        Path to the tool executable if found, None otherwise.
    """
    return shutil.which(tool_name)


def verify_tool_dependencies(tool_type: ToolType) -> None:
    """
    Verify that the required external tool is installed.
    
    Args:
        tool_type (ToolType): The type of tool to check for.
    
    Raises:
        RuntimeError: If the tool is not found in PATH.
    """
    tool_map = {
        ToolType.FOLDSEEK: "foldseek",
        ToolType.TMALIGN: "TMalign",
        ToolType.USALIGN: "USalign"
    }
    
    tool_name = tool_map[tool_type]
    tool_path = check_external_tool(tool_name)
    
    if tool_path is None:
        raise RuntimeError(
            f"Required tool '{tool_name}' not found in PATH.\n"
            f"Please install {tool_name} and ensure it's accessible.\n"
            f"Installation instructions: [add relevant links]"
        )
    
    print(f"Found {tool_name} at: {tool_path}")