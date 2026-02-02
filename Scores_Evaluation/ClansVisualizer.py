import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


class ClansVisualizer:
    """Handles visualization operations for CLANS data."""
    
    @staticmethod
    def generate_scatter_plot(
        data_x: pd.Series, 
        data_y: pd.Series, 
        x_label: str = "X-axis", 
        y_label: str = "Y-axis", 
        title: str = "Scatter Plot", 
        save_path: str | None = None
    ):
        """
        Generates a scatter plot comparing two sets of values.

        Args:
            data_x (pd.Series or list-like): Values for the x-axis.
            data_y (pd.Series or list-like): Values for the y-axis.
            x_label (str, optional): Label for x-axis. Defaults to "X-axis".
            y_label (str, optional): Label for y-axis. Defaults to "Y-axis".
            title (str, optional): Plot title. Defaults to "Scatter Plot".
            save_path (str, optional): Path to save the figure. If None, plot is shown. Defaults to None.

        Returns:
            None
        """
        data_x = pd.Series(data_x)
        data_y = pd.Series(data_y)
        
        if len(data_x) != len(data_y):
            raise ValueError("data_x and data_y must have the same length.")

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

    @staticmethod
    def generate_bar_plot(
        data: pd.DataFrame,
        x_label: str = "X-axis",
        y_label: str = "Y-axis",
        title: str = "Bar Plot",
        save_path: str | None = None,
        is_grouped: bool = False,
        bar_labels: list[str] | None = None
    ):
        """
        Generates a bar plot from a DataFrame.

        Args:
            data (pd.DataFrame): DataFrame where index is x-axis and columns are y-axis values.
            x_label (str, optional): Label for x-axis. Defaults to "X-axis".
            y_label (str, optional): Label for y-axis. Defaults to "Y-axis".
            title (str, optional): Plot title. Defaults to "Bar Plot".
            save_path (str, optional): Path to save the figure. If None, plot is shown. Defaults to None.
            is_grouped (bool): If True, creates a grouped bar plot. Defaults to False.
            bar_labels (list[str], optional): Labels for the bars in a grouped plot. Defaults to None.
        """
        fig, ax = plt.subplots(figsize=(10, 6))
        
        if is_grouped:
            if bar_labels:
                data[bar_labels].plot(kind='bar', ax=ax)
            else:
                data.plot(kind='bar', ax=ax)
        else:
            data.plot(kind='bar', ax=ax, legend=False)

        ax.set_title(title, fontsize=16)
        ax.set_ylabel(y_label, fontsize=12)
        ax.set_xlabel(x_label, fontsize=12)
        ax.tick_params(axis='x', rotation=45)
        if is_grouped:
            ax.legend(title='Category')

        plt.tight_layout()

        if save_path:
            plt.savefig(save_path, dpi=300)
            plt.close()
        else:
            plt.show()

    @staticmethod
    def generate_heatmap(
        pivot_table: pd.DataFrame,
        title: str = "Heatmap",
        x_label: str = "X-axis",
        y_label: str = "Y-axis",
        save_path: str | None = None,
        annot: bool = True,
        cmap: str = 'viridis',
        fmt: str = ".2f"
    ):
        """
        Generates a heatmap from a pivot table.

        Args:
            pivot_table (pd.DataFrame): Pivot table to visualize.
            title (str, optional): Plot title. Defaults to "Heatmap".
            x_label (str, optional): Label for x-axis. Defaults to "X-axis".
            y_label (str, optional): Label for y-axis. Defaults to "Y-axis".
            save_path (str, optional): Path to save the figure. If None, plot is shown. Defaults to None.
            annot (bool): If True, write the data value in each cell. Defaults to True.
            cmap (str): The mapping from data values to color space. Defaults to 'viridis'.
            fmt (str): String formatting code to use when adding annotations. Defaults to ".2f".
        """
        plt.figure(figsize=(12, 8))
        sns.heatmap(pivot_table, annot=annot, cmap=cmap, fmt=fmt)
        plt.title(title, fontsize=16) 
        plt.xlabel(x_label, fontsize=12)
        plt.ylabel(y_label, fontsize=12)
        plt.tight_layout()
        if save_path:
            plt.savefig(save_path, dpi=300)
            plt.close()
        else:
            plt.show()
