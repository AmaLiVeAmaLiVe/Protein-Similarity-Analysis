"""
msa.py
Progressive Multiple Sequence Alignment (MSA) engine.
Progressively aligns sequences to maintain a uniform coordinate register.
"""

from typing import List
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from src.aligner import get_configured_aligner


def progressive_msa(records: List[SeqRecord]) -> List[SeqRecord]:
    """
    Performs progressive alignment using the first sequence as the seed,
    progressively expanding the alignment grid so all sequences have identical lengths.

    Args:
        records: List of Bio.SeqRecord objects.

    Returns:
        List of Bio.SeqRecord objects with equal sequence lengths (including gaps).
    """
    if not records:
        return []

    aligner = get_configured_aligner()

    # Initialize aligned pool with the first record
    seed_rec = records[0]
    aligned_seqs = [str(seed_rec.seq)]
    record_ids = [seed_rec.id]

    # Progressively align each incoming record against the current consensus/seed
    for target_rec in records[1:]:
        tgt_id = target_rec.id
        tgt_raw_seq = str(target_rec.seq)

        # Align reference seed (ungapped) to target sequence (ungapped)
        ref_ungapped = aligned_seqs[0].replace("-", "")
        aln = aligner.align(ref_ungapped, tgt_raw_seq)[0]
        aln_ref_str = str(aln[0])
        aln_tgt_str = str(aln[1])

        # Map new gaps into all existing aligned sequences
        new_aligned_seqs = [[] for _ in aligned_seqs]
        new_target_seq = []

        ref_char_idx = 0
        for r_char, t_char in zip(aln_ref_str, aln_tgt_str):
            new_target_seq.append(t_char)

            if r_char != "-":
                # Advance existing aligned sequences by one ungapped character
                for s_idx, seq_str in enumerate(aligned_seqs):
                    # Find corresponding position in existing sequence
                    char_count = 0
                    pos_found = None
                    for pos, ch in enumerate(seq_str):
                        if ch != "-":
                            if char_count == ref_char_idx:
                                pos_found = pos
                                break
                            char_count += 1

                    if pos_found is not None:
                        new_aligned_seqs[s_idx].append(seq_str[pos_found])
                    else:
                        new_aligned_seqs[s_idx].append("-")

                ref_char_idx += 1
            else:
                # A gap was opened in the reference: insert gaps across all existing sequences
                for s_idx in range(len(aligned_seqs)):
                    new_aligned_seqs[s_idx].append("-")

        # Update aligned pool
        aligned_seqs = ["".join(chars) for chars in new_aligned_seqs]
        aligned_seqs.append("".join(new_target_seq))
        record_ids.append(tgt_id)

    # Pad any trailing length mismatches to ensure a perfectly rectangular matrix
    max_len = max(len(s) for s in aligned_seqs)
    padded_seqs = [s.ljust(max_len, "-") for s in aligned_seqs]

    # Package back into Bio.SeqRecord objects
    msa_result = [
        SeqRecord(Seq(s), id=rec_id, description="")
        for s, rec_id in zip(padded_seqs, record_ids)
    ]

    return msa_result


