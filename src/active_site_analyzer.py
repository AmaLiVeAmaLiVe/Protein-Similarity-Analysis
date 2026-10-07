"""
active_site_analyzer.py
Enhanced universal analyzer that maps functional residues across homologous proteins.
Includes local window search (tolerance for alignment drift/indels) and 
physicochemical property evaluation.
"""

from typing import List, Dict, Tuple, Optional
import pandas as pd
from Bio.SeqRecord import SeqRecord
from src.aligner import get_configured_aligner

# Groupings of amino acids by biochemical properties
AMINO_ACID_GROUPS = {
    "Acidic": {"D", "E"},
    "Basic": {"K", "R", "H"},
    "Polar_Neutral": {"S", "T", "N", "Q"},
    "Hydrophobic_Aliphatic": {"A", "V", "L", "I", "M"},
    "Hydrophobic_Aromatic": {"F", "Y", "W"},
    "Special": {"P", "G", "C"},
}

def get_aa_group(aa: str) -> str:
    for group, members in AMINO_ACID_GROUPS.items():
        if aa in members:
            return group
    return "Unknown"


def map_reference_residue_enhanced(
    ref_seq: str, 
    target_seq: str, 
    ref_pos_1indexed: int, 
    expected_aa: str, 
    aligner,
    window_tolerance: int = 2
) -> Dict[str, any]:
    """
    Maps a 1-indexed position from ref_seq to target_seq.
    If the exact aligned position doesn't match expected_aa, it checks
    within +/- window_tolerance to detect alignment drift around gaps.
    """
    alignment = aligner.align(ref_seq, target_seq)[0]
    aligned_ref = str(alignment[0])
    aligned_tgt = str(alignment[1])

    current_ref_pos = 0
    current_tgt_pos = 0
    target_aligned_idx = None
    target_pos_at_site = None
    target_aa_at_site = "-"

    for idx, (r_char, t_char) in enumerate(zip(aligned_ref, aligned_tgt)):
        if r_char != "-":
            current_ref_pos += 1
        if t_char != "-":
            current_tgt_pos += 1

        if current_ref_pos == ref_pos_1indexed and r_char != "-":
            target_aligned_idx = idx
            target_aa_at_site = t_char
            target_pos_at_site = current_tgt_pos if t_char != "-" else None
            break

    # Exact match check
    exact_match = (target_aa_at_site == expected_aa)
    shifted_match_found = False
    shifted_pos = None

    # If exact match failed and tolerance is enabled, scan local sequence window
    if not exact_match and target_pos_at_site is not None:
        start_scan = max(1, target_pos_at_site - window_tolerance)
        end_scan = min(len(target_seq), target_pos_at_site + window_tolerance)
        
        for candidate_pos in range(start_scan, end_scan + 1):
            if candidate_pos != target_pos_at_site and target_seq[candidate_pos - 1] == expected_aa:
                shifted_match_found = True
                shifted_pos = candidate_pos
                break

    # Determine status
    if exact_match:
        status = "Conserved"
    elif shifted_match_found:
        status = f"Conserved (Shifted to pos {shifted_pos})"
    else:
        # Check if mutation preserves biochemical group
        exp_group = get_aa_group(expected_aa)
        tgt_group = get_aa_group(target_aa_at_site)
        if exp_group == tgt_group:
            status = f"Conservative ({exp_group})"
        else:
            status = "Non-conserved"

    return {
        "residue": f"{target_aa_at_site}{target_pos_at_site if target_pos_at_site else '-'}",
        "exact_conserved": exact_match,
        "shifted_pos": shifted_pos,
        "status": status,
    }


def inspect_active_sites(
    records: List[SeqRecord],
    ref_id: str,
    sites_config: Dict[str, Dict[str, any]],
    window_tolerance: int = 2
) -> pd.DataFrame:
    """
    Enhanced active-site analyzer supporting tolerance for alignment drift.
    """
    aligner = get_configured_aligner()
    ref_rec = next((r for r in records if r.id == ref_id), None)
    if not ref_rec:
        raise ValueError(f"Reference '{ref_id}' not found.")

    ref_seq = str(ref_rec.seq)

    # Validate reference positions
    for site_name, cfg in sites_config.items():
        pos = cfg["pos"]
        exp_aa = cfg["expected_aa"]
        actual_aa = ref_seq[pos - 1]
        if actual_aa != exp_aa:
            raise ValueError(
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
                row[f"{site_name}_Status"] = "Reference"
                row[f"{site_name}_Conserved"] = True
                continue

            res_info = map_reference_residue_enhanced(
                ref_seq=ref_seq,
                target_seq=tgt_seq,
                ref_pos_1indexed=pos,
                expected_aa=exp_aa,
                aligner=aligner,
                window_tolerance=window_tolerance
            )

            row[f"{site_name}_Residue"] = res_info["residue"]
            row[f"{site_name}_Status"] = res_info["status"]
            # Conserved if exact match OR shifted nearby match found
            row[f"{site_name}_Conserved"] = res_info["exact_conserved"] or (res_info["shifted_pos"] is not None)

        rows.append(row)

    return pd.DataFrame(rows)
