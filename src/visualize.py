"""Generalizes heatmaps and dendrograms from similarity matrices."""

from pathlib import Path
import matplotlib.pyplot as plt 
import pandas as pd
import seaborn as sns
from scipy.cluster import hierarchy


def plot_heatmap(df_identity: pd.DataFrame, output_path: str | Path):
    """Plots and saves an annotated pairwise identity heatmap"""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fig_size = max(8, len(df_identity)*0.7)
    plt.figure(figsize=(fig_size, fig_size*0.8))

    sns.heatmap(
        df_identity,
        annot=True,
        fmt=".1f",
        cmap="viridis",
        linewidths=0.5,
        cbar_kws={"label": "% Sequence Identity"},
    )
    plt.title("Pairwise Protein Sequence Identity (BLOSUM62)", fontsize=14, pad=15)
    plt.xticks(rotation=45, ha="right")
    plt.yticks(rotation=0)
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def plot_dendrogram(linkage_matrix, labels: list[str], output_path: str | Path):
    """Plots and saves a horizontal hierarchical clustering dendrogram"""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    fig_height = max(5, len(labels)*0.4)
    plt.figure(figsize=(10, fig_height))

    hierarchy.dendrogram(
        linkage_matrix,
        labels=labels,
        orientation="left",
        leaf_font_size=10,
    )
    plt.title("Hierarchical Clustering (Sequence Divergence)", fontsize=14, pad=15)
    plt.xlabel("Distance Score (100 - % Identity)", fontsize=11)
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()

