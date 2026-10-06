"""Computes global pairwise alignments and similarity/distance metrics."""


from typing import List, Tuple
import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.SeqRecord import SeqRecord


def get_configured_aligner() -> Align.PairwiseAligner:
    """Initializes a PairWiseAligner with BLOSUM62 and standard affine gap costs."""
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5

    return aligner


def calculate_identity(aligned_seq1: str, aligned_seq2: str, min_len: int) -> float:
    """
    Calculates exact percentage identity between two aligned sequences.
    Positions where both contain gaps are ignores; mismatches/indels are penalized.
    """
    matches = sum(
        res1 == res2
        for res1, res2 in zip(aligned_seq1, aligned_seq2)
        if res1 != "-" and res2 != "-"
    )

    return (matches/min_len) * 100.00


def compute_similiraity_matrix(records: List[SeqRecord]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Builds symmetric Percent Identity and Distance matrices for all sequence pairs.
    Args:
        record: list of parced SeqRecord objects.
    Returns:
        Tuple of (df_identity, df_distance) as pandas.DataFrames.
    """
    aligner = get_configured_aligner()
    n = len(records)
    names = [rec.id for rec in records]
    identity_matrix = np.zeros((n, n), dtype=float)

    for i in range(n):
        identity_matrix[i, i] = 100.0     # Perferct self-identity
        seq1 = str(records[i].seq)

        for j in range(i+1, n):
            seq2 = str(records[j].seq)
            alignment = aligner.align(seq1, seq2)[0]

            # alignment[0] and alignment[1] contain string representations of aligned chains
            pid = calculate_identity(
                alignment[0], alignment[1], min_len=min(len(seq1), len(seq2))
                )
            identity_matrix[i, j] = pid
            identity_matrix[j, i] = pid

    distance_matrix = 100.0 - identity_matrix

    df_identity = pd.DataFrame(identity_matrix, index=names, columns=names)
    df_distance = pd.DataFrame(distance_matrix, index=names, columns=names)

    return df_identity, df_distance
