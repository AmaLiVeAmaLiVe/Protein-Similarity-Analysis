"""
Phylogenetic tree construction and analysis module.
Constructs UPGMA and Neighbourhood-Joining (NJ) trees from distance matrices,
computes branch lengths, exports Newick format, and generates tree figures.
"""


from pathlib import Path
from typing import Optional, Union
import matplotlib.pyplot as plt
import pandas as pd
from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor


def dataframe_to_biopython_dm(df_distance: pd.DataFrame) -> DistanceMatrix:
    """
    Converts a pandas symmetric N x N distance DataFrame into a 
    BioPython DistanceMatrix object (lower-triangular matrix) 
    """
    names = list(df_distance.index)
    matrix = []

    for i, name in enumerate(names):
        # DistanceMatrix expects lower triangle: row i with columns 0..i
        row = [float(df_distance.loc[name, names[j]]) for j in range(i+1)]
        matrix.append(row)

    return DistanceMatrix(names, matrix)


def build_phylogenetic_tree(df_distance: pd.DataFrame, method: str = "nj") -> Tree:
    """
    Constructs a phylogenetic tree from a distance matrix
    Args:
        df_distance: N x N evolutionary distance matrix (e.g., 100 - %ID)
        method: 'nj' (Neighbourhood-Joining, default), or 'upgma'
    Returns:
        Bio.Phylo Tree object with calculated branch lengths 
    """
    dm = dataframe_to_biopython_dm(df_distance)
    constructor = DistanceTreeConstructor()

    if method.lower() == "nj":
        tree = constructor.nj(dm)
    elif method.lower() in ("upgma", "ugpma"):
        tree = constructor.upgma(dm)
    else:
        raise ValueError(f"Unknown tree method '{method}'. Choose 'nj' or 'upgma'.")

    # Clean internal node names for clean Newick rendering
    for clade in tree.find_clades():
        if clade.name and clade.name.startswith("Inner"):
            clade.name = None

    return tree


def export_newick(tree: Tree, output_path: Union[str, Path]):
    """Exports the phylogenetic tree to a standard Newick (.nwk) file"""
    path = Path(output_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    Phylo.write(tree, str(path), "newick")


def plot_phylogenetic_tree(
        tree: Tree,
        title: str = "Neighbourhood-Joining Phylogenetic Tree",
        output_path: Optional[Union[str, Path]] = None
):
    """
    Renders the phylogenetic tree with explicit branch lengths
    """
    fig, ax = plt.subplots(figsize=(9, 4.5))
    Phylo.draw(
        tree,
        axes=ax,
        do_show=False,
        show_confidence=False,
        branch_labels=lambda c: f"{c.branch_length:.2f}" if c.branch_length else ""
    )
    ax.set_title(title, fontsize=13, weight="bold", pad=12)
    ax.set_xlabel("Evolutionary Distance (Substitutions / 100 aa)")
    plt.tight_layout()

    if output_path:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(out, dpi=300)

    return fig
    

