from __future__ import annotations

from pathlib import Path

from .ensembl import SPECIES
from .sequence import clean_dna


TUBB5_GENE_ID = "ENSMUSG00000001525"
TUBB5_TRANSCRIPT_ID = "ENSMUST00000001566"
# The bundled sequence is the UCSC/GENCODE VM snapshot carrying transcript version .10.
# Online mode resolves the current Ensembl version of the same stable transcript.
TUBB5_FIXTURE_VERSION = "10"
TUBB5_CURRENT_CANONICAL_VERSION_AT_BUILD = "11"
TUBB5_CHROMOSOME = "17"
TUBB5_STRAND = -1
TUBB5_CDNA_LENGTH = 2649
TUBB5_CDS_START_CDNA0 = 253
TUBB5_INSERTION_CDNA0 = 1585
TUBB5_STOP_CODON = "TAA"
TUBB5_INSERTION_BOUNDARY0 = 36_145_876
TUBB5_STOP_START0 = 36_145_873
TUBB5_STOP_END0 = 36_145_876
TUBB5_PROTEIN_LENGTH_AA = 444


def load_tubb5_cdna() -> str:
    path = Path(__file__).resolve().parents[1] / "data" / "tubb5_ucsc_mrna.fa"
    lines = [line.strip() for line in path.read_text(encoding="utf-8").splitlines()]
    sequence = clean_dna("".join(line for line in lines if line and not line.startswith(">")))
    if len(sequence) != TUBB5_CDNA_LENGTH:
        raise ValueError(
            f"Bundled Tubb5 fixture has length {len(sequence)}, expected {TUBB5_CDNA_LENGTH}"
        )
    return sequence


TUBB5_SPECIES = SPECIES["mouse"]
