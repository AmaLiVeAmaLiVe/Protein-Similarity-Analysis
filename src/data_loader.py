"""
Parses and validates protein FASTA files.
"""


from pathlib import Path
from typing import List
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# Standard 20 amino acid single-letter codes plus common IUPAC characters
VALID_AMINO_ACIDS = set("ACDEFGHIKLMNPQRSTVWYBXZU*")


def clean_header(record_id: str) -> str:
    """Extract a concise identifier from typical UniProt or NCBI headers"""
    # Example: 'sp|P69905|HBA_HUMAN' -> 'HBA_HUMAN'
    parts = record_id.split("|")

    return parts[-1] if len(parts) > 1 else record_id


def load_sequences(filepath: str | Path) -> List[SeqRecord]:
    """
    Loads, validates, and normalizes sequences from a FASTA file.
    Args:
        filepath: Path to the input FASTA file
    Returns:
        List of cleaned SeqRecord objects
    """
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"FASTA file not found at: {path.resolve()}")

    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise ValueError(f"No FASTA sequences found in: {path}")

    cleaned_records = []
    seend_ids = set()

    for rec in records:
        short_id = clean_header(rec.id)
        if short_id in seend_ids:
            raise ValueError(f"Duplicate sequence identifier detected: {short_id}")

        seend_ids.add(short_id)

        seq_str = str(rec.seq).upper()
        invalid_chars = set(seq_str) - VALID_AMINO_ACIDS
        if invalid_chars:
            raise ValueError(
                f"Sequence '{short_id}' contains invalid amino acid characters: {invalid_chars}"
            )

        rec.id = short_id
        rec.seq = rec.seq.upper()
        cleaned_records.append(rec)

    return cleaned_records 
