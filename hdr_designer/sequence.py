from __future__ import annotations

from collections.abc import Iterable

_DNA_COMPLEMENT = str.maketrans("ACGTRYMKBDHVNacgtrymkbdhvn", "TGCAYRKMVHDBNtgcayrkmvhd bn".replace(" ", ""))

# Standard genetic code. Kept local to avoid a heavyweight dependency.
CODON_TABLE: dict[str, str] = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def clean_dna(sequence: str) -> str:
    """Return uppercase DNA with whitespace removed and validate IUPAC symbols."""
    cleaned = "".join(sequence.split()).upper()
    allowed = set("ACGTRYMKBDHVN")
    invalid = sorted(set(cleaned) - allowed)
    if invalid:
        raise ValueError(f"Invalid DNA symbol(s): {', '.join(invalid)}")
    return cleaned


def reverse_complement(sequence: str) -> str:
    return clean_dna(sequence).translate(_DNA_COMPLEMENT)[::-1]


def gc_percent(sequence: str) -> float:
    sequence = clean_dna(sequence)
    if not sequence:
        return 0.0
    gc = sequence.count("G") + sequence.count("C")
    return round(100.0 * gc / len(sequence), 1)


def translate(sequence: str, *, stop_at_stop: bool = False) -> str:
    sequence = clean_dna(sequence)
    if len(sequence) % 3:
        raise ValueError("Coding sequence length is not divisible by three")
    peptide: list[str] = []
    for index in range(0, len(sequence), 3):
        codon = sequence[index:index + 3]
        amino_acid = CODON_TABLE.get(codon, "X")
        if amino_acid == "*" and stop_at_stop:
            break
        peptide.append(amino_acid)
    return "".join(peptide)


def wrap_sequence(sequence: str, width: int = 60) -> str:
    sequence = clean_dna(sequence)
    return "\n".join(sequence[index:index + width] for index in range(0, len(sequence), width))


def fasta_record(name: str, sequence: str, width: int = 60) -> str:
    safe_name = " ".join(name.strip().split())
    return f">{safe_name}\n{wrap_sequence(sequence, width)}\n"


def find_motif_positions(sequence: str, motifs: Iterable[str]) -> dict[str, list[int]]:
    sequence = clean_dna(sequence)
    result: dict[str, list[int]] = {}
    for raw_motif in motifs:
        motif = clean_dna(raw_motif)
        positions: list[int] = []
        start = 0
        while True:
            found = sequence.find(motif, start)
            if found < 0:
                break
            positions.append(found)
            start = found + 1
        result[motif] = positions
    return result
