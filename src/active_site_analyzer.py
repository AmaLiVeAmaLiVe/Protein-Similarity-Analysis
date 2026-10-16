"""
Universal analyzer that maps arbitrary functional residues across homologous
sequences and evaluates both strict conservation and physicochemical shifts.
"""


from typing import List, Dict, Tuple, Optional
import pandas as pd
from Bio.SeqRecord import SeqRecord
from src.aligner import get_configured_aligner


# Standard IUPAC amino acid physicochemical classifications
AA_PROPERTIES = {
    # Negatively Charged / Acidic
    "D": "Acidic (-)", "E": "Acidic (-)",
    # Positively Charged / Basic
    "K": "Basic (+)", "R": "Basic (+)", "H": "Basic (+)",
    # Polar Uncharged / Neutral
    "S": "Polar Neutral", "T": "Polar Neutral", "N": "Polar Neutral", "Q": "Polar Neutral",
    # Hydrophobic Aliphatic
    "A": "Hydrophobic Aliphatic", "V": "Hydrophobic Aliphatic", 
    "L": "Hydrophobic Aliphatic", "I": "Hydrophobic Aliphatic", "M": "Hydrophobic Aliphatic",
    # Hydrophobic Aromatic
    "F": "Aromatic", "Y": "Aromatic", "W": "Aromatic",
    # Special Conformation
    "P": "Conformational (Pro)", "G": "Conformational (Gly)", "C": "Crosslinking (Cys)"
}


def get_property(aa: str) -> str:
    """Returns the biochemical class of an amino acid."""
    return AA_PROPERTIES.get(aa.upper(), "Unknown")


def map_reference_residue(
    ref_seq: str, 
    target_seq: str, 
    ref_pos_1indexed: int, 
    expected_aa: str, 
    aligner=None,
    window_tolerance: int = 2
) -> Dict[str, any]:
    """
    Maps a 1-indexed reference coordinate to target_seq.
    Uses adaptive Needleman-Wunsch / Smith-Waterman routing.
    """
    from src.aligner import get_configured_aligner, select_alignment_mode

    # Auto-select mode based on length ratio
    mode = select_alignment_mode(ref_seq, target_seq)
    active_aligner = get_configured_aligner(mode=mode)
    
    alignment = active_aligner.align(ref_seq, target_seq)[0]
    
    # In Biopython PairwiseAlignment, coordinates track exact start/end slices
    aligned_ref = str(alignment[0])
    aligned_tgt = str(alignment[1])
    ref_start = alignment.coordinates[0][0]
    tgt_start = alignment.coordinates[1][0]

    current_ref_pos = ref_start
    current_tgt_pos = tgt_start
    target_pos = None
    target_aa = "-"

    for r_char, t_char in zip(aligned_ref, aligned_tgt):
        if r_char != "-":
            current_ref_pos += 1
        if t_char != "-":
            current_tgt_pos += 1

        if current_ref_pos == ref_pos_1indexed and r_char != "-":
            target_aa = t_char
            target_pos = current_tgt_pos if t_char != "-" else None
            break

    # Tolerance fallback & property classification
    exact_match = (target_aa == expected_aa)
    shifted_pos = None

    if not exact_match and target_pos is not None:
        start_scan = max(1, target_pos - window_tolerance)
        end_scan = min(len(target_seq), target_pos + window_tolerance)
        for candidate_pos in range(start_scan, end_scan + 1):
            if candidate_pos != target_pos and target_seq[candidate_pos - 1] == expected_aa:
                shifted_pos = candidate_pos
                break

    ref_prop = get_property(expected_aa)
    tgt_prop = get_property(target_aa)

    if exact_match:
        interpretation = "Conserved Identical"
        is_conserved = True
    elif shifted_pos is not None:
        interpretation = f"Conserved (Register shifted to pos {shifted_pos})"
        is_conserved = True
    elif target_aa == "-":
        interpretation = "Deleted / Gapped"
        is_conserved = False
    elif ref_prop == tgt_prop:
        interpretation = f"Conservative Swap ({ref_prop})"
        is_conserved = True
    else:
        interpretation = f"Diverged ({ref_prop} -> {tgt_prop})"
        is_conserved = False

    residue_label = f"{target_aa}{target_pos if target_pos is not None else '-'}"

    return {
        "residue": residue_label,
        "is_conserved": is_conserved,
        "interpretation": interpretation,
        "alignment_mode": mode
    }


def inspect_active_sites(
    records: List[SeqRecord],
    ref_id: str,
    sites_config: Dict[str, Dict[str, any]],
    window_tolerance: int = 2
) -> pd.DataFrame:
    """
    Universal functional residue inspector.
    """
    aligner = get_configured_aligner()
    ref_rec = next((r for r in records if r.id == ref_id), None)
    if not ref_rec:
        raise ValueError(f"Reference '{ref_id}' not found among loaded sequences.")

    ref_seq = str(ref_rec.seq)

    # Validate reference positions against expected amino acids
    for site_name, cfg in sites_config.items():
        pos = cfg["pos"]
        exp_aa = cfg["expected_aa"]
        actual_aa = ref_seq[pos - 1]
        if actual_aa != exp_aa:
            raise ValueError(
                f"Reference validation error for '{site_name}': "
                f"Position {pos} in '{ref_id}' is '{actual_aa}', expected '{exp_aa}'."
            )

    rows = []
    for rec in records:
        row = {"Protein": rec.id}
        tgt_seq = str(rec.seq)

        for site_name, cfg in sites_config.items():
            pos = cfg["pos"]
            exp_aa = cfg["expected_aa"]

            if rec.id == ref_id:
                row[f"{site_name}_Residue"] = f"{exp_aa}{pos}"
                row[f"{site_name}_Effect"] = "Reference"
                row[f"{site_name}_Conserved"] = True
                continue

            res_info = map_reference_residue(
                ref_seq=ref_seq,
                target_seq=tgt_seq,
                ref_pos_1indexed=pos,
                expected_aa=exp_aa,
                aligner=aligner,
                window_tolerance=window_tolerance
            )

            row[f"{site_name}_Residue"] = res_info["residue"]
            row[f"{site_name}_Effect"] = res_info["interpretation"]
            row[f"{site_name}_Conserved"] = res_info["is_conserved"]

        rows.append(row)

    return pd.DataFrame(rows)
