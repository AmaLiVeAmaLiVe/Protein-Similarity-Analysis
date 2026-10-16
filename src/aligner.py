"""
Configures pairwise alignment algorithms with deterministic auto-routing
between Needleman-Wunsch (global) and Smith-Waterman (local).
"""

from typing import Tuple, List, Optional
import numpy as np
import pandas as pd
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.SeqRecord import SeqRecord


def get_configured_aligner(mode: str = "global") -> Align.PairwiseAligner:
    """
    Initializes a Biopython PairwiseAligner with standard BLOSUM62 parameters.
    mode can be 'global' (Needleman-Wunsch) or 'local' (Smith-Waterman).
    """
    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10.0
    aligner.extend_gap_score = -0.5
    return aligner


def select_alignment_mode(
    seq1: str, 
    seq2: str, 
    length_ratio_threshold: float = 0.75
) -> str:
    """
    Deterministically selects alignment strategy:
    - If length ratio >= threshold: 'global' (Needleman-Wunsch)
    - If length ratio < threshold:  'local' (Smith-Waterman)
    """
    l1, l2 = len(seq1), len(seq2)
    ratio = min(l1, l2) / max(l1, l2)
    return "global" if ratio >= length_ratio_threshold else "local"


def align_sequences_adaptive(
    seq1: str, 
    seq2: str,
    forced_mode: Optional[str] = None
) -> Tuple[str, str, str]:
    """
    Aligns two sequences using the auto-selected or specified algorithm.
    Returns: (aligned_seq1, aligned_seq2, chosen_mode)
    """
    mode = forced_mode if forced_mode else select_alignment_mode(seq1, seq2)
    aligner = get_configured_aligner(mode=mode)
    alignment = aligner.align(seq1, seq2)[0]
    return str(alignment[0]), str(alignment[1]), mode


def calculate_identity(aligned_seq1: str, aligned_seq2: str, min_len: int) -> float:
    """
    Calculates exact percentage identity between two aligned sequences.
    Positions where both contain gaps are ignored; mismatches/indels are penalized.
    """
    matches = sum(
        res1 == res2
        for res1, res2 in zip(aligned_seq1, aligned_seq2)
        if res1 != "-" and res2 != "-"
    )
    return (matches / min_len) * 100.00


def compute_similarity_matrix(
    records: List[SeqRecord],
    adaptive: bool = True
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Computes NxN percent identity and evolutionary distance matrices.
    If adaptive=True, automatically chooses global vs local alignment per pair.
    """
    n = len(records)
    ids = [r.id for r in records]
    id_matrix = np.zeros((n, n), dtype=float)
    dist_matrix = np.zeros((n, n), dtype=float)

    for i in range(n):
        id_matrix[i, i] = 100.0
        dist_matrix[i, i] = 0.0
        for j in range(i + 1, n):
            s1 = str(records[i].seq)
            s2 = str(records[j].seq)

            aln_s1, aln_s2, _ = align_sequences_adaptive(s1, s2) if adaptive else (
                *align_sequences_adaptive(s1, s2, forced_mode="global")[:2], "global"
            )

            # Calculate identity over aligned non-gap matches
            matches = sum(1 for a, b in zip(aln_s1, aln_s2) if a == b and a != "-")
            min_len = min(len(s1), len(s2))
            pct_identity = (matches / min_len) * 100.0

            id_matrix[i, j] = pct_identity
            id_matrix[j, i] = pct_identity
            
            dist = max(0.0, 100.0 - pct_identity)
            dist_matrix[i, j] = dist
            dist_matrix[j, i] = dist

    df_id = pd.DataFrame(id_matrix, index=ids, columns=ids)
    df_dist = pd.DataFrame(dist_matrix, index=ids, columns=ids)
    return df_id, df_dist
