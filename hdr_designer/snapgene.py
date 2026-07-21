from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from pathlib import Path
from typing import Any

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


class SnapGeneError(RuntimeError):
    """Raised when a SnapGene file cannot be read or fails basic validation."""


@dataclass(frozen=True)
class SnapGeneDocument:
    path: Path
    record: SeqRecord
    sequence: str
    sha256: str
    topology: str
    features: tuple[dict[str, Any], ...]

    @property
    def length(self) -> int:
        return len(self.sequence)


def _feature_label(feature: Any) -> str:
    for key in ("label", "gene", "product", "note"):
        values = feature.qualifiers.get(key)
        if values:
            value = values[0] if isinstance(values, list) else values
            return str(value)
    return feature.type


def _serialize_feature(feature: Any) -> dict[str, Any]:
    parts = [
        {
            "start0": int(part.start),
            "end0": int(part.end),
            "strand": int(part.strand or 0),
        }
        for part in feature.location.parts
    ]
    return {
        "type": str(feature.type),
        "label": _feature_label(feature),
        "start0": min(part["start0"] for part in parts),
        "end0": max(part["end0"] for part in parts),
        "strand": int(feature.location.strand or 0),
        "parts": parts,
    }


def read_snapgene(path: str | Path) -> SnapGeneDocument:
    """Read a SnapGene .dna file through Biopython's binary parser."""
    file_path = Path(path)
    if not file_path.is_file():
        raise SnapGeneError(f"SnapGene file not found: {file_path}")
    try:
        record = SeqIO.read(str(file_path), "snapgene")
    except Exception as exc:  # Biopython raises several parser-specific exceptions.
        raise SnapGeneError(f"Could not parse SnapGene file {file_path.name}: {exc}") from exc

    sequence = str(record.seq).upper()
    invalid = sorted(set(sequence) - set("ACGTN"))
    if invalid:
        raise SnapGeneError(
            f"SnapGene sequence contains unsupported symbols: {', '.join(invalid)}"
        )
    if not sequence:
        raise SnapGeneError("SnapGene file contains no DNA sequence")

    topology = str(record.annotations.get("topology", "linear")).lower()
    return SnapGeneDocument(
        path=file_path,
        record=record,
        sequence=sequence,
        sha256=sha256(file_path.read_bytes()).hexdigest(),
        topology=topology,
        features=tuple(_serialize_feature(feature) for feature in record.features),
    )
