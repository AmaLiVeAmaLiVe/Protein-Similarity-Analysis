"""Performs hierarchical linkage clustering on protein distance matrices."""

import numpy as np
import pandas as pd
from scipy.cluster import hierarchy
from scipy.spatial.distance import squareform


def compute_hierarchical_linkage(
        df_distance: pd.DataFrame, method: str="average"
        ) -> np.ndarray:
    """
    Calculates hierarchical linkage from a symmetric distance DataFrame.
    Args:
        df_distance: Square symmetric distance DataFrame
        method: Linkage algortihm ("average" for UPGMA, "complete", or "ward")
    Returns:
        Linkage matrics as NumPy array for dendrogram plotting
    """
    matrix = df_distance.to_numpy(copy=True)

    # Ensure numerical symmetry on a main diagonal
    np.fill_diagonal(matrix, 0.0)

    # Convert square matrix to condensed form (upper triangle) required by SciPy
    condensed_dist = squareform(matrix)
    linkage_matrix = hierarchy.linkage(condensed_dist, method=method)

    return linkage_matrix
