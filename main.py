#!/usr/bin/env python3
"""
main.py
CLI entry point orchestrating the protein similarity pipeline:
1. Loads and validates FASTA sequences.
2. Computes pairwise alignment matrices (identity and distance).
3. Computes hierarchical clustering and exports visualizations.
4. (Optional) Analyzes active site conservation if a JSON config is provided.
"""

import argparse
import json
import sys
from pathlib import Path

from src.data_loader import load_sequences
from src.aligner import compute_similarity_matrix
from src.clustering import compute_hierarchical_linkage
from src.visualize import plot_heatmap, plot_dendrogram
from src.active_site_analyzer import inspect_active_sites


def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Protein Pairwise Similarity, Divergence, and Active Site Analysis Pipeline"
    )
    parser.add_argument(
        "--input", "-i",
        type=Path,
        required=True,
        help="Path to the input protein FASTA file.",
    )
    parser.add_argument(
        "--output-dir", "-o",
        type=Path,
        default=Path("results"),
        help="Root folder for output tables and figures (default: results).",
    )
    parser.add_argument(
        "--sites-config", "-c",
        type=Path,
        default=None,
        help="Optional path to JSON configuration defining reference protein and active sites.",
    )
    return parser.parse_args()


def run_active_site_inspection(records, config_path: Path, output_table: Path) -> None:
    """Helper function to load config, inspect residues, and save the result table."""
    if not config_path.is_file():
        print(f"Warning: Sites config file not found at {config_path}. Skipping active site analysis.")
        return

    with open(config_path, "r", encoding="utf-8") as f:
        try:
            config_data = json.load(f)
        except json.JSONDecodeError as exc:
            print(f"Error parsing JSON in {config_path}: {exc}. Skipping active site analysis.")
            return

    ref_id = config_data.get("reference_id")
    sites = config_data.get("sites")

    if not ref_id or not isinstance(sites, dict):
        print("Error: Sites config must contain 'reference_id' (str) and 'sites' (dict). Skipping.")
        return

    print(f"Running active site conservation analysis against reference '{ref_id}'...")
    try:
        df_sites = inspect_active_sites(
            records=records,
            ref_id=ref_id,
            sites_config=sites
        )
        df_sites.to_csv(output_table, index=False)
        print(f"Saved active site conservation table to: {output_table}")
        print("\nActive Site Conservation Summary:")
        print(df_sites.to_string(index=False))
        print()
    except Exception as exc:
        print(f"Error during active site inspection: {exc}")


def main():
    args = parse_arguments()

    if not args.input.is_file():
        sys.exit(f"Error: Input sequence file does not exist: {args.input}")

    # Set up directory layout
    tables_dir = args.output_dir / "tables"
    figures_dir = args.output_dir / "figures"
    tables_dir.mkdir(parents=True, exist_ok=True)
    figures_dir.mkdir(parents=True, exist_ok=True)

    # 1. Data Ingestion
    print(f"Loading sequences from: {args.input}")
    records = load_sequences(args.input)
    print(f"Loaded {len(records)} valid protein sequences.")

    # 2. Pairwise Alignment & Similarity Matrices
    print("Computing global pairwise alignments (BLOSUM62)...")
    df_identity, df_distance = compute_similarity_matrix(records)

    identity_csv = tables_dir / "identity_matrix.csv"
    distance_csv = tables_dir / "distance_matrix.csv"
    df_identity.to_csv(identity_csv)
    df_distance.to_csv(distance_csv)
    print(f"Saved similarity matrices to: {tables_dir}/")

    # 3. Hierarchical Clustering and Visualization
    print("Generating clustering models and visualizations...")
    linkage = compute_hierarchical_linkage(df_distance, method="average")

    heatmap_path = figures_dir / "identity_heatmap.png"
    dendrogram_path = figures_dir / "tree_dendrogram.png"

    plot_heatmap(df_identity, heatmap_path)
    plot_dendrogram(linkage, list(df_identity.index), dendrogram_path)
    print(f"Visualizations saved to: {figures_dir}/")

    # 4. Active Site Analysis (only runs if the flag was provided)
    if args.sites_config is not None:
        sites_csv = tables_dir / "active_sites_conservation.csv"
        run_active_site_inspection(records, args.sites_config, sites_csv)

    print("Pipeline execution complete.")


if __name__ == "__main__":
    main()
    