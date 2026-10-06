#!/usr/bin/env python3
"""CLI entrypoint to run the protein similarity analysis pipeline"""


import argparse
from pathlib import Path
from src.data_loader import load_sequences
from src.aligner import compute_similiraity_matrix
from src.clustering import compute_hierarchical_linkage
from src.visualize import plot_heatmap, plot_dendrogram


def main():
    parser = argparse.ArgumentParser(
        description="Protein Pairwise Similarity and Divergence Analysis Pipeline"
    )
    parser.add_argument(
        "--input",
        type=str,
        default="data/raw/sequences.fasta",
        help="Path to the input protein FASTA file.",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="results",
        help="Root folder for output tables and figures.",
    )
    args = parser.parse_args()

    input_path = Path(args.input)
    output_base = Path(args.output_dir)
    tables_dir = output_base / "tables"
    figures_dir = output_base / "figures"

    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    print(f"Loading sequences from: {input_path}")
    records = load_sequences(input_path)
    print(f"Loaded {len(records)} valid protein sequences.")

    print("Computing global pairwise alignments (BLOSUM62)...")
    df_identity, df_distance = compute_similiraity_matrix(records)

    # Save matrix data
    identity_csv = tables_dir / "identity_matrix.csv"
    distance_csv = tables_dir / "distance_matrix.csv"
    df_identity.to_csv(identity_csv)
    df_distance.to_csv(distance_csv)
    print(f"Saved matrices to: {tables_dir}/")

    print("Generating clustering models and visualizations...")
    linkage = compute_hierarchical_linkage(df_distance, method="average")

    heatmap_path = figures_dir / "identity_heatmap.png"
    dendrogram_path = figures_dir / "tree_dendrogram.png"

    plot_heatmap(df_identity, heatmap_path)
    plot_dendrogram(linkage, list(df_identity.index), dendrogram_path)

    print(f"Visualizations successfully saved to: {figures_dir}/")
    print("Pipeline execution complete.")


if __name__ == "__main__":
    main()
