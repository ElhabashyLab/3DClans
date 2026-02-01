import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class ClansVisualizer:
    """Handles visualization operations for CLANS data."""
    
    @staticmethod
    def generate_scatter_plot(
        data_x: pd.Series, 
        data_y: pd.Series, 
        x_label: str | None = None, 
        y_label: str | None = None, 
        title: str | None = None, 
        save_path: str | None = None
    ):
        """
        Generates a scatter plot comparing two sets of values.

        Args:
            data_x (pd.Series or list-like): Values for the x-axis.
            data_y (pd.Series or list-like): Values for the y-axis.
            x_label (str, optional): Label for x-axis. Defaults to None.
            y_label (str, optional): Label for y-axis. Defaults to None.
            title (str, optional): Plot title. Defaults to None.
            save_path (str, optional): Path to save the figure. If None, plot is shown. Defaults to None.

        Returns:
            None
        """
        data_x = pd.Series(data_x)
        data_y = pd.Series(data_y)
        
        if len(data_x) != len(data_y):
            raise ValueError("data_x and data_y must have the same length.")

        # Default labels if not provided
        if x_label is None:
            x_label = "X-axis"
        if y_label is None:
            y_label = "Y-axis"
        if title is None:
            title = f"{y_label} vs {x_label}"

        plt.figure(figsize=(6, 6))
        sns.scatterplot(x=data_x, y=data_y, alpha=0.7, edgecolor=None)
        
        # Plot identity line y=x
        min_val = min(data_x.min(), data_y.min())
        max_val = max(data_x.max(), data_y.max())
        plt.plot([min_val, max_val], [min_val, max_val], linestyle='--', color='gray', alpha=0.5)

        plt.xlabel(x_label)
        plt.ylabel(y_label)
        plt.title(title)
        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300)
            plt.close()
        else:
            plt.show()
