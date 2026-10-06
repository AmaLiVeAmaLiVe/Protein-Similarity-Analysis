"""
active_site_analyzer.py
Universal analyzer that maps arbitrary functional residues from any reference
protein across homologous target sequences using global alignment.
"""

from typing import List, Dict, Optional, Tuple
import pandas as pd
from Bio import Align
from Bio.SeqRecord import SeqRecord
from src.aligner import get_configured_aligner


def map_residue(
    ref_seq: str, target_seq: str, ref_pos: int, aligner: Align.PairwiseAligner
) -> Tuple[str, Optional[int]]:
    """
    Finds which residue in target_seq aligns with position ref_pos (1-indexed) in ref_seq.

    Returns:
        (target_residue, target_1indexed_pos) or ('-', None) if aligned to a gap.
    """
    # Force global alignment
    alignments = aligner.align(ref_seq, target_seq)
    if not alignments:
        return "-", None

    alignment = alignments[0]
    
    # Biopython alignment string format: line 0 is target/ref, line 2 is the other
    # Using format(alignment, "fasta") guarantees standardized string extraction
    fasta_blocks = format(alignment, "fasta").split(">")
    # Block 1 is ref, Block 2 is target
    seq_ref_aligned = "".join(fasta_blocks[1].split("\n")[1:]).strip().upper()
    seq_target_aligned = "".join(fasta_blocks[2].split("\n")[1:]).strip().upper()

    ref_idx = 0      # 1-indexed count in original reference
    target_idx = 0   # 1-indexed count in original target

    for r_aa, t_aa in zip(seq_ref_aligned, seq_target_aligned):
        if r_aa != "-":
            ref_idx += 1
        if t_aa != "-":
            target_idx += 1

        if ref_idx == ref_pos:
            if r_aa == "-":
                # Position is inside an insertion relative to reference
                continue
            if t_aa == "-":
                return "-", None
            return t_aa, target_idx

    return "-", None


def inspect_active_sites(
    records: List[SeqRecord],
    ref_id: str,
    sites_config: Dict[str, Dict[str, any]],
) -> pd.DataFrame:
    """
    Universal active-site conservation analyzer.
    """
    aligner = get_configured_aligner()
    ref_record = next((r for r in records if r.id == ref_id), None)
    if not ref_record:
        raise ValueError(f"Reference '{ref_id}' not found in loaded sequences.")

    ref_seq = str(ref_record.seq).upper()
    results = []

    # Validate reference positions against configuration
    for site_name, cfg in sites_config.items():
        pos = cfg["pos"]
        exp_aa = cfg["expected_aa"]
        if pos < 1 or pos > len(ref_seq):
            raise ValueError(f"Position {pos} out of range for {ref_id} (len: {len(ref_seq)})")
        actual_ref_aa = ref_seq[pos - 1]
        if actual_ref_aa != exp_aa:
            raise ValueError(
                f"Reference validation error for '{site_name}': "
                f"Position {pos} in '{ref_id}' is '{actual_ref_aa}', but expected '{exp_aa}'."
            )

    for rec in records:
        row = {"Protein": rec.id}
        target_seq = str(rec.seq).upper()

        if rec.id == ref_id:
            for site_name, cfg in sites_config.items():
                pos = cfg["pos"]
                aa = cfg["expected_aa"]
                row[f"{site_name}_Residue"] = f"{aa}{pos}"
                row[f"{site_name}_Conserved"] = True
            results.append(row)
            continue

        for site_name, cfg in sites_config.items():
            ref_pos = cfg["pos"]
            exp_aa = cfg["expected_aa"]
            target_aa, target_pos = map_residue(ref_seq, target_seq, ref_pos, aligner)

            pos_str = str(target_pos) if target_pos is not None else "-"
            row[f"{site_name}_Residue"] = f"{target_aa}{pos_str}"
            row[f"{site_name}_Conserved"] = (target_aa == exp_aa)

        results.append(row)

    return pd.DataFrame(results)
