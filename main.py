#!/usr/bin/env python3
"""
main.py
Comprehensive CLI entry point orchestrating:
1. FASTA Loading & Pairwise Needleman-Wunsch Alignment.
2. UPGMA Evolutionary Dendrogram & Identity Heatmap.
3. Progressive Multiple Sequence Alignment (MSA).
4. Shannon Entropy Conservation Profiling.
5. Active Site Inspection.
"""

import argparse
from pathlib import Path
import sys

from src.data_loader import load_sequences
from src.aligner import compute_similarity_matrix
from src.clustering import compute_hierarchical_linkage
from src.visualize import plot_heatmap, plot_dendrogram
from src.active_site_analyzer import inspect_active_sites
from src.msa import progressive_msa
from src.entropy import compute_msa_conservation, plot_conservation_profile


def main():
    parser = argparse.ArgumentParser(description="Full Protein Similarity & Conservation Pipeline")
    parser.add_argument("--input", "-i", type=Path, default=Path("data/raw/sequences.fasta"))
    parser.add_argument("--output-dir", "-o", type=Path, default=Path("results"))
    parser.add_argument("--sites-config", "-c", type=Path, default=Path("config/sites/globins.json"))
    args = parser.parse_args()

    tables_dir = args.output_dir / "tables"
    figures_dir = args.output_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    # 1. Load Data
    print(f"[1/5] Loading sequences from: {args.input}")
    records = load_sequences(args.input)
    print(f"      Loaded {len(records)} valid protein sequences.")

    # 2. Pairwise Similarity & UPGMA
    print("[2/5] Computing pairwise alignments & distance tree...")
    df_id, df_dist = compute_similarity_matrix(records)
    df_id.to_csv(tables_dir / "identity_matrix.csv")
    df_dist.to_csv(tables_dir / "distance_matrix.csv")

    plot_heatmap(df_id, figures_dir / "identity_heatmap.png")
    linkage = compute_hierarchical_linkage(df_dist, method="average")
    plot_dendrogram(linkage, list(df_id.index), figures_dir / "tree_dendrogram.png")

    # 3. Progressive Multiple Sequence Alignment (MSA)
    print("[3/5] Constructing Progressive Multiple Sequence Alignment (MSA)...")
    msa_records = progressive_msa(records)

    # 4. Shannon Entropy & Conservation Profiling
    print("[4/5] Computing per-position Shannon Entropy profile...")
    df_entropy = compute_msa_conservation(msa_records)
    df_entropy.to_csv(tables_dir / "shannon_entropy_profile.csv", index=False)

    # In globin MSA, Distal and Proximal correspond to alignment columns ~65 and ~94
    highlights = [(65, "Distal His"), (94, "Proximal His")]
    plot_conservation_profile(df_entropy, highlights, figures_dir / "conservation_profile.png")

    # 5. Active Site Analysis
    print("[5/5] Performing active site tolerance analysis...")
    if args.sites_config.is_file():
        import json
        with open(args.sites_config, "r", encoding="utf-8") as f:
            cfg = json.load(f)
        df_sites = inspect_active_sites(records, cfg["reference_id"], cfg["sites"])
        df_sites.to_csv(tables_dir / "active_sites_conservation.csv", index=False)
        print("\nActive Site Summary:")
        print(df_sites.to_string(index=False))

    print(f"\nAll artifacts generated successfully in: {args.output_dir.resolve()}")


if __name__ == "__main__":
    main()




    