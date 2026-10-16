"""
Phylogenetic tree construction and analysis module.
Constructs UPGMA and Neighbourhood-Joining (NJ) trees from distance matrices,
computes branch lengths, exports Newick format, and generates tree figures.
"""


from pathlib import Path
from typing import Optional, Union, Tuple, List, Dict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from Bio import Phylo
from Bio.Phylo.BaseTree import Tree
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
import plotly.graph_objects as go


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


def _get_radial_coords(
    clade, 
    x_coords: Dict, 
    y_coords: Dict, 
    radial_coords: Dict,
    total_taxa: int
):
    """
    Transforms rectangular tree coordinates (x=root distance, y=leaf index)
    into circular polar coordinates (theta in radians, r=radius).
    """
    x = x_coords[clade]
    y = y_coords[clade]
    
    # Map vertical range 0..total_taxa to angle 0..2*pi
    theta = (2.0 * np.pi * y) / max(1, total_taxa)
    r = x  # Evolutionary distance from root is radius
    
    x_cart = r * np.cos(theta)
    y_cart = r * np.sin(theta)
    
    radial_coords[clade] = (x_cart, y_cart, r, np.degrees(theta))


def render_interactive_tree(
    tree: Phylo.BaseTree.Tree,
    layout: str = "circular",
    title: str = "Interactive Phylogenetic Tree",
    width: int = 750,
    height: int = 700
) -> go.Figure:
    """
    Generates an interactive Plotly figure from a Bio.Phylo Tree object.
    
    Args:
        tree: Bio.Phylo Tree object.
        layout: 'circular' (radial cladogram) or 'rectangular'.
        title: Chart title.
        width: Canvas width in pixels.
        height: Canvas height in pixels.
        
    Returns:
        go.Figure with hoverable tooltips, pan, and zoom.
    """
    # 1. Compute tree coordinates
    clades = list(tree.find_clades(order="level"))
    terminals = tree.get_terminals()
    total_taxa = len(terminals)

    # Calculate depth (distance from root)
    depths = tree.depths()
    
    # Assign y-coordinates to terminals
    y_coords = {}
    for idx, term in enumerate(terminals):
        y_coords[term] = idx + 1
        
    # Assign y-coordinates to internal nodes (average of children)
    for clade in tree.get_nonterminals(order="postorder"):
        child_ys = [y_coords[c] for c in clade]
        y_coords[clade] = sum(child_ys) / len(child_ys)

    radial_coords = {}
    for clade in clades:
        _get_radial_coords(clade, depths, y_coords, radial_coords, total_taxa)

    # 2. Build lines (branches)
    edge_x = []
    edge_y = []
    
    for parent in tree.get_nonterminals():
        px, py, pr, pth = radial_coords[parent] if layout == "circular" else (depths[parent], y_coords[parent], 0, 0)
        for child in parent:
            cx, cy, cr, cth = radial_coords[child] if layout == "circular" else (depths[child], y_coords[child], 0, 0)
            
            edge_x.extend([px, cx, None])
            edge_y.extend([py, cy, None])

    # 3. Build Nodes & Leaf Labels
    node_x = []
    node_y = []
    hover_texts = []
    labels_text = []
    node_colors = []
    node_sizes = []

    for clade in clades:
        is_leaf = clade.is_terminal()
        if layout == "circular":
            x, y, r, deg = radial_coords[clade]
        else:
            x, y = depths[clade], y_coords[clade]
            r = x

        node_x.append(x)
        node_y.append(y)

        br_len = f"{clade.branch_length:.2f}" if clade.branch_length is not None else "0.00"
        name = clade.name if clade.name else "Internal Node"
        
        hover = (
            f"<b>{name}</b><br>"
            f"Branch Length: {br_len}<br>"
            f"Distance to Root: {r:.2f}<br>"
            f"Type: {'Terminal Leaf' if is_leaf else 'Ancestral Split'}"
        )
        hover_texts.append(hover)
        
        labels_text.append(f"  {clade.name}" if is_leaf else "")
        node_colors.append("#e74c3c" if is_leaf else "#3498db")
        node_sizes.append(10 if is_leaf else 6)

    # 4. Construct Plotly Figure
    fig = go.Figure()

    # Add branch lines
    fig.add_trace(go.Scatter(
        x=edge_x,
        y=edge_y,
        mode="lines",
        line=dict(color="#7f8c8d", width=2),
        hoverinfo="none",
        showlegend=False
    ))

    # Add nodes + text labels
    fig.add_trace(go.Scatter(
        x=node_x,
        y=node_y,
        mode="markers+text",
        marker=dict(size=node_sizes, color=node_colors, line=dict(width=1, color="#2c3e50")),
        text=labels_text,
        textposition="top right",
        textfont=dict(size=12, family="Arial", color="#2c3e50"),
        hoverinfo="text",
        hovertext=hover_texts,
        showlegend=False
    ))

    axis_visible = (layout != "circular")
    fig.update_layout(
        title=dict(text=title, font=dict(size=15, family="Arial")),
        width=width,
        height=height,
        template="plotly_white",
        xaxis=dict(showgrid=axis_visible, zeroline=False, showticklabels=axis_visible, title="Substitutions / 100 aa" if axis_visible else ""),
        yaxis=dict(showgrid=axis_visible, zeroline=False, showticklabels=axis_visible, title="Taxa" if axis_visible else ""),
        margin=dict(l=40, r=60, b=40, t=50)
    )

    return fig

