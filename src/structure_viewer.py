"""
structure_viewer.py
Fetches 3D coordinates from RCSB PDB and builds a py3Dmol visualizer
highlighting the heme pocket and active site histidines.
"""

from pathlib import Path
import urllib.request
import py3Dmol


def fetch_pdb_file(pdb_id: str, save_dir: str | Path = "data/structures") -> str:
    """
    Downloads .pdb file if not already present, returns the PDB text content.
    """
    save_path = Path(save_dir) / f"{pdb_id.lower()}.pdb"
    save_path.parent.mkdir(parents=True, exist_ok=True)

    if not save_path.exists():
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        req = urllib.request.Request(
            url,
            headers={"User-Agent": "BioinformaticsProject/1.0"}
        )
        with urllib.request.urlopen(req) as resp:
            content = resp.read().decode("utf-8")
        with open(save_path, "w", encoding="utf-8") as f:
            f.write(content)
    else:
        with open(save_path, "r", encoding="utf-8") as f:
            content = f.read()

    return content


def render_active_sites_3d(
    pdb_id: str = "1A3N",
    chain: str = "A",
    distal_pos: int = 58,
    proximal_pos: int = 87
) -> py3Dmol.view:
    """
    Renders an interactive 3D view of the globin monomer, highlighting
    the heme cofactor and the dual catalytic histidines.
    """
    pdb_text = fetch_pdb_file(pdb_id)
    if not pdb_text or len(pdb_text) < 100:
        raise ValueError(f"Failed to load valid PDB data for '{pdb_id}'.")

    # Create viewer with explicit canvas dimensions
    view = py3Dmol.view(width=750, height=500)
    view.addModel(pdb_text, "pdb")

    # 1. Background style: Dim secondary chains (B, C, D)
    view.setStyle(
        {"chain": ["B", "C", "D"]},
        {"cartoon": {"color": "#d3d3d3", "opacity": 0.25}}
    )

    # 2. Focus chain (Chain A): Soft cyan cartoon backbone
    view.setStyle(
        {"chain": chain},
        {"cartoon": {"color": "#4682b4", "opacity": 0.9}}
    )

    # 3. Heme group (HEM in Chain A or general HEM)
    # Shown as thick orange sticks with central iron sphere
    view.addStyle(
        {"resn": "HEM"},
        {"stick": {"colorscheme": "orangeCarbon", "radius": 0.25}}
    )
    view.addStyle(
        {"resn": "HEM", "atom": "FE"},
        {"sphere": {"color": "#ff3300", "radius": 1.2}}
    )

    # 4. Catalytic Histidines (Distal & Proximal)
    # Shown as bright green sticks with all non-hydrogen atoms visible
    view.addStyle(
        {"chain": chain, "resi": [distal_pos, proximal_pos]},
        {"stick": {"colorscheme": "greenCarbon", "radius": 0.3}}
    )

    # 5. Position labels
    view.addLabel(
        f"Distal His{distal_pos}",
        {"chain": chain, "resi": distal_pos},
        {"backgroundColor": "#222222", "fontColor": "#00ff00", "fontSize": 12, "backgroundOpacity": 0.8}
    )
    view.addLabel(
        f"Proximal His{proximal_pos}",
        {"chain": chain, "resi": proximal_pos},
        {"backgroundColor": "#222222", "fontColor": "#00ff00", "fontSize": 12, "backgroundOpacity": 0.8}
    )

    # 6. Camera focus directly on the active site pocket
    view.zoomTo({"chain": chain, "resi": [distal_pos, proximal_pos], "resn": "HEM"})
    view.render()

    return view




