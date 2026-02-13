import pandas as pd
import numpy as np


class DataNormalizer:
    """Handles all data normalization operations."""
    
    @staticmethod
    def normalize(df: pd.DataFrame, column_name: str, normalization_type: str) -> pd.Series:
        """
        Normalizes a specified column in the DataFrame using the given normalization type.

        Args:
            df (pd.DataFrame): The DataFrame containing the column to normalize.
            column_name (str): The name of the column to normalize.
            normalization_type (str): The type of normalization to apply. Supported types: ("min-max", "z-score", "-log10").

        Returns:
            pd.Series: A Series with the normalized column.
        """
        if column_name not in df.columns:
            raise ValueError(f"Column '{column_name}' not found in DataFrame.")
        
        col = df[column_name].astype(float)
        
        normalizers = {
            "min-max": DataNormalizer._normalize_min_max,
            "z-score": DataNormalizer._normalize_z_score,
            "-log10": DataNormalizer._normalize_log10
        }
        
        if normalization_type not in normalizers:
            raise ValueError(
                f"Unsupported normalization_type '{normalization_type}'. "
                "Supported: 'min-max', 'z-score', '-log10'."
            )
        
        return normalizers[normalization_type](col)
    
    @staticmethod
    def _normalize_min_max(col: pd.Series) -> pd.Series:
        """
        Applies min-max normalization to a pandas Series.
        Scales all values to the range [0, 1]. If the column is constant
        (max == min), all values are set to 0.

        Args:
            col (pd.Series): The input numeric column to normalize.

        Returns:
            pd.Series: A Series with min-max normalized values in [0, 1].
        """
        min_val = col.min()
        max_val = col.max()
        if max_val == min_val:
            return pd.Series(0.0, index=col.index)
        return (col - min_val) / (max_val - min_val)

    @staticmethod
    def _normalize_z_score(col: pd.Series) -> pd.Series:
        """
        Applies z-score normalization (standardization) to a pandas Series.
        Each value is transformed as (x - mean) / std. If the column has zero
        variance (std == 0), all values are set to 0.

        Args:
            col (pd.Series): The input numeric column to normalize.

        Returns:
            pd.Series: A Series with z-score normalized values.
        """
        std = col.std()
        if std == 0:
            return pd.Series(0.0, index=col.index)
        return (col - col.mean()) / std

    @staticmethod
    def _normalize_log10(col: pd.Series) -> pd.Series:
        """
        Applies -log10 transformation to a pandas Series making small values more interpretable.
        Raises a ValueError if the column contains zero or negative values.

        Args:
            col (pd.Series): The input numeric column to transform. Must be > 0.

        Returns:
            pd.Series: A Series where each value is transformed as -log10(x).

        Raises:
            ValueError: If any value in the column is <= 0.
        """
        if (col <= 0).any():
            raise ValueError("Cannot apply -log10 to zero or negative values.")
        return pd.Series(-np.log10(col.values), index=col.index)
