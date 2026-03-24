import enum

"""This module defines the ToolType enum class,
 which represents the different types of similarity tools that can be used in the 3DClans pipeline."""
class ToolType(enum.Enum):
    FOLDSEEK = "foldseek"
    USALIGN = "USalign"
