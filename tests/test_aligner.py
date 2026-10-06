"""
test_aligner.py
Unit tests verifying alignment logic and matrix outputs.
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from src.aligner import compute_similarity_matrix, calculate_identity


def test_calculate_identity_exact():
    seq1 = "HEAGAWGHEE"
    seq2 = "HEAGAWGHEE"
    assert calculate_identity(seq1, seq2, len(seq1)) == 100.0


def test_calculate_identity_with_mismatch():
    seq1 = "HEAGAWGHEE"
    seq2 = "HEAGAWGHAE"  # 1 mismatch at index 8
    assert calculate_identity(seq1, seq2, len(seq1)) == 90.0


def test_matrix_dimensions_and_diagonal():
    records = [
        SeqRecord(Seq("MVLSPADKTN"), id="P1"),
        SeqRecord(Seq("MVHLTPEEKS"), id="P2"),
    ]
    df_id, df_dist = compute_similarity_matrix(records)

    assert df_id.shape == (2, 2)
    assert df_id.loc["P1", "P1"] == 100.0
    assert df_id.loc["P2", "P2"] == 100.0
    assert df_dist.loc["P1", "P1"] == 0.0
    assert df_id.loc["P1", "P2"] == df_id.loc["P2", "P1"]
    