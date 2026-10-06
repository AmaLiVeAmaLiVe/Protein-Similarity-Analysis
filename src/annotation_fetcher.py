"""
annotation_fetcher.py
Fetches experimentally verified functional annotations (active sites,
binding sites, metal coordination) directly from the UniProtKB REST API.
"""

import json
from pathlib import Path
from typing import Dict, Any
import urllib.request
import urllib.error


def fetch_uniprot_features(uniprot_id: str, output_path: str | Path | None = None) -> Dict[str, Any]:
    """
    Downloads curated active sites and binding sites for a UniProt ID.
    
    Args:
        uniprot_id: UniProt accession (e.g., 'P69905') or entry name ('HBA_HUMAN').
        output_path: Optional path to save the generated JSON config.
        
    Returns:
        Structured dictionary matching our pipeline's sites configuration.
    """
    url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}.json"
    req = urllib.request.Request(
        url,
        headers={"User-Agent": "ProteinSimilarityPipeline/1.0 (academic; student-project)"}
    )

    try:
        with urllib.request.urlopen(req) as response:
            data = json.loads(response.read().decode("utf-8"))
    except urllib.error.HTTPError as exc:
        raise RuntimeError(f"Failed to fetch UniProt data for '{uniprot_id}': HTTP {exc.code}")

    sequence = data["sequence"]["value"]
    entry_id = data.get("uniProtkbId", uniprot_id)
    features = data.get("features", [])

    # Filter for biologically active sites and ligand binding residues
    target_feature_types = {"Active site", "Binding site"}
    extracted_sites = {}

    for idx, feat in enumerate(features):
        ft_type = feat.get("type")
        if ft_type in target_feature_types:
            start_pos = feat["location"]["start"]["value"]
            end_pos = feat["location"]["end"]["value"]

            # Only track single-residue functional sites
            if start_pos == end_pos:
                res_aa = sequence[start_pos - 1]
                desc = feat.get("description", ft_type).replace(" ", "_")
                # Clean label for key
                site_key = f"{ft_type.replace(' ', '_')}_{start_pos}_{desc[:25]}"
                
                extracted_sites[site_key] = {
                    "pos": start_pos,
                    "expected_aa": res_aa,
                    "description": feat.get("description", "")
                }

    config_payload = {
        "reference_id": entry_id,
        "accession": data.get("primaryAccession", uniprot_id),
        "source": "UniProtKB/Swiss-Prot",
        "sites": extracted_sites
    }

    if output_path:
        out_p = Path(output_path)
        out_p.parent.mkdir(parents=True, exist_ok=True)
        with open(out_p, "w", encoding="utf-8") as f:
            json.dump(config_payload, f, indent=2)
        print(f"Curated annotations saved to: {out_p.resolve()}")

    return config_payload


if __name__ == "__main__":
    import sys
    accession = sys.argv[1] if len(sys.argv) > 1 else "P69905"
    out_file = sys.argv[2] if len(sys.argv) > 2 else "config/sites/uniprot_P69905.json"
    print(f"Fetching verified functional annotations for UniProt ID: {accession}...")
    fetch_uniprot_features(accession, out_file)
    