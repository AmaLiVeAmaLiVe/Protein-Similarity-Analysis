from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from src.active_site_analyzer import inspect_active_sites


def test_physicochemical_divergence_reporting():
    # Synthetic test sequences:
    # Ref has Glutamate (E, Acidic) at pos 2
    # Target 1 has Aspartate (D, Acidic) -> Conservative Swap
    # Target 2 has Valine (V, Hydrophobic) -> Diverged
    records = [
        SeqRecord(Seq("MEAAAA"), id="REF"),
        SeqRecord(Seq("MDAAAA"), id="CONSERVATIVE_TARGET"),
        SeqRecord(Seq("MVAAAA"), id="DIVERGED_TARGET"),
    ]

    sites_cfg = {
        "Test_Site": {"pos": 2, "expected_aa": "E"}
    }

    df = inspect_active_sites(records, "REF", sites_cfg)

    # Check Conservative Swap detection
    row_cons = df[df["Protein"] == "CONSERVATIVE_TARGET"].iloc[0]
    assert row_cons["Test_Site_Conserved"] == True
    assert "Conservative Swap" in row_cons["Test_Site_Effect"]

    # Check Diverged detection
    row_div = df[df["Protein"] == "DIVERGED_TARGET"].iloc[0]
    assert row_div["Test_Site_Conserved"] == False
    assert "Diverged (Acidic (-) -> Hydrophobic Aliphatic)" in row_div["Test_Site_Effect"]
    