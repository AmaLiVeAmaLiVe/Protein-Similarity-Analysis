"""
Calculates per-column Shannon entropy and plots conservation profiles across MSAs. 
"""


from collections import Counter
import math
from pathlib import Path
from typing import List, Tuple
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord


def calculate_column_entropy(column: List[str]) -> float:
    """
    Computes Shannon entropy (base 2) for an amino acid alignment column.
    Gaps are excluded from frequency counts.
    """
    residues = [aa for aa in column if aa != "-"]
    if not residues:
        return 0.0

    total = len(residues)
    counts = Counter(residues)
    entropy = 0.0

    for count in counts.values():
        p = count/total
        entropy -= p*math.log2(p)

    return entropy


def compute_msa_conservation(msa_records: List[SeqRecord]) -> pd.DataFrame:
    """Computes Shannon entropy and normalized conservation (0.0 to 1.0) for every column"""
    num_cols = len(msa_records[0].seq)
    max_possible_entropy = math.log2(20)     # Max diversity across 20 amino acids (~4.32 bits)

    rows = []
    for col_idx in range(num_cols):
        col_residues = [str(r.seq)[col_idx] for r in msa_records]
        H = calculate_column_entropy(col_residues)
        normalized_conservation = max(0.0, 1.0-(H / max_possible_entropy))

        # Track consensus residue
        counts = Counter([aa for aa in col_residues if aa != "-"])
        dominant_aa = counts.most_common(1)[0][0] if counts else "-"

        rows.append({
            "Alignment_Position": col_idx+1,
            "Dominant_Residue": dominant_aa,
            "Shannon_Entropy": round(H, 4),
            "Conservation_Score": round(normalized_conservation, 4)
        })

    return pd.DataFrame(rows)


def plot_conservation_profile(
        df_conservation: pd.DataFrame,
        highlight_positions: List[Tuple[int, str]],
        ouput_path: str | Path
):
    """
    Plots conservation score across sequence positions with active site highlighted
    """
    path = Path(ouput_path)
    path.parent.mkdir(parents=True, exist_ok=True)

    plt.figure(figsize=(14, 4.5))
    plt.plot(
        df_conservation["Alignment_Position"],
        df_conservation["Conservation_Score"],
        color="#1f77b4",
        linewidth=1.5,
        label="Per-position Conservation"
    )

    # Highligh active sites
    for pos, label in highlight_positions:
        score = df_conservation.loc[
            df_conservation["Alignment_Position"] == pos, "Conservation_Score"
        ].values

        if len(score) > 0:
            plt.scatter(pos, score[0], color="#d62728", s=80, zorder=5)
            plt.annotate(
                f"{label} (Pos {pos})",
                xy=(pos, score[0]),
                xytext=(pos + 2, score[0] - 0.15),
                arrowprops=dict(facecolor="#d62728", shrink=0.05, width=1, headwidth=6),
                fontsize=9,
                fontweight="bold"
            )

    plt.title("Per-Residue Conservation Profile Across Globin MSA", fontsize=13, pad=12)
    plt.xlabel("Alignment Column Position", fontsize=11)
    plt.ylabel("Conservation Index (1 - Normalized H)", fontsize=11)
    plt.ylim(-0.05, 1.1)
    plt.grid(True, linestyle="--", alpha=0.6)
    plt.legend(loc="lower right")
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()

