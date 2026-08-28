from __future__ import annotations

import csv
import hashlib
import io
import re
from typing import Any
import zipfile

from .exports import (
    assembled_plasmid_genbank,
    genotyping_primers_csv,
    locus_context_genbank,
)
from .generate_oligos_from_guides import guide_oligos_csv
from .models import DesignResult, HomologyArm
from .synthesis_qc import (
    TWIST_MAX_HOMOPOLYMER_NT,
    TWIST_ORDERING_RULESET,
    homopolymer_findings,
)


ORDERING_GENBANK_DATE = "01-JAN-1980"


class OrderingError(ValueError):
    """Raised when a design is not ready to be packaged for ordering."""


def _arm_ordering_findings(arm: HomologyArm) -> list[dict[str, Any]]:
    return [
        {"arm": arm.name, **finding}
        for finding in homopolymer_findings(arm.final_gene_oriented_sequence)
    ]


def twist_ordering_qc(result: DesignResult) -> dict[str, Any]:
    """Check the final mutation-containing homology arms for Twist ordering."""
    findings = [
        *_arm_ordering_findings(result.five_prime_arm),
        *_arm_ordering_findings(result.three_prime_arm),
    ]
    adjustments = [
        {"arm": arm.name, **arm.boundary_adjustment}
        for arm in (result.five_prime_arm, result.three_prime_arm)
        if arm.boundary_adjustment
    ]
    return {
        "status": "ERROR" if findings else "PASS",
        "ruleset": TWIST_ORDERING_RULESET,
        "max_homopolymer_nt": TWIST_MAX_HOMOPOLYMER_NT,
        "findings": findings,
        "boundary_adjustments": adjustments,
        "detail": (
            f"Found {len(findings)} homopolymer run(s) longer than "
            f"{TWIST_MAX_HOMOPOLYMER_NT} nt in the final homology arms."
            if findings
            else (
                "Neither final homology arm contains a homopolymer longer than "
                f"{TWIST_MAX_HOMOPOLYMER_NT} nt."
            )
        ),
        "vendor_screening_note": (
            "This local check applies the configured homopolymer limit only. "
            "Twist's ordering portal must still screen the complete submitted sequences."
        ),
    }


def design_identity(result: DesignResult) -> str:
    identity_fields = (
        result.assembly,
        result.transcript_id,
        result.terminus,
        result.backbone_addgene_id,
        result.guides[0].spacer if result.guides else "no-guide",
        result.five_prime_arm.final_gene_oriented_sequence,
        result.three_prime_arm.final_gene_oriented_sequence,
        str(result.donor_payload.get("payload_sequence_5to3", "")),
    )
    return hashlib.sha256("\x1f".join(identity_fields).encode("utf-8")).hexdigest()[:10]


def _safe_stem(result: DesignResult) -> str:
    tag_name = (
        result.backbone_name
        if result.backbone_addgene_id == "custom"
        else str(result.donor_payload.get("tag_name", "tag"))
    )
    raw = (
        f"{result.gene_symbol}_{result.terminus.lower().replace('-', '_')}_"
        f"{tag_name}_"
        f"{design_identity(result)}"
    )
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", raw).strip("._") or "hdr_design"


def ordering_package_filename(result: DesignResult) -> str:
    return f"{_safe_stem(result)}_ordering_package.zip"


def _require_ordering_ready(result: DesignResult) -> dict[str, Any]:
    if not result.sequence_complete:
        raise OrderingError("The design is not sequence-complete, so no ordering package can be released.")
    if not result.locus_contexts:
        raise OrderingError("Annotated wild-type and edited locus records are unavailable.")
    qc = twist_ordering_qc(result)
    if qc["status"] != "PASS":
        details = "; ".join(
            f"{item['arm']} bases {item['interval_1based']} "
            f"({item['base']} x {item['length_nt']})"
            for item in qc["findings"]
        )
        raise OrderingError(
            "Twist homopolymer validation failed: final homology arms may not contain "
            f"runs longer than {TWIST_MAX_HOMOPOLYMER_NT} nt. {details}"
        )
    return qc


def twist_sequences_csv(result: DesignResult) -> str:
    """Return the two final SapI-flanked homology-arm fragments for Twist upload."""
    qc = _require_ordering_ready(result)
    stem = _safe_stem(result)
    rows = (
        ("UHA", result.cloning_fragments.get("uha_synthesis_fragment_5to3")),
        ("DHA", result.cloning_fragments.get("dha_synthesis_fragment_5to3")),
    )
    if any(not sequence for _, sequence in rows):
        raise OrderingError("One or both final homology-arm synthesis fragments are unavailable.")

    fields = [
        "Name",
        "Sequence",
        "Length",
        "Sequence Type",
        "Gene",
        "Terminus",
        "SHA-256",
        "Internal QC Status",
        "Internal QC Ruleset",
        "Requested Arm Length",
        "Final Arm Length",
        "Boundary Adjustment",
        "Twist Portal Screening",
    ]
    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    arms = (result.five_prime_arm, result.three_prime_arm)
    for (arm_label, sequence_value), arm in zip(rows, arms):
        sequence = str(sequence_value).upper()
        writer.writerow(
            {
                "Name": f"{stem}_{arm_label}",
                "Sequence": sequence,
                "Length": len(sequence),
                "Sequence Type": f"final_{arm_label.lower()}_synthesis_fragment",
                "Gene": result.gene_symbol,
                "Terminus": result.terminus,
                "SHA-256": hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                "Internal QC Status": qc["status"],
                "Internal QC Ruleset": qc["ruleset"],
                "Requested Arm Length": arm.requested_length or arm.length,
                "Final Arm Length": arm.length,
                "Boundary Adjustment": arm.correction_note,
                "Twist Portal Screening": "REQUIRED",
            }
        )
    return output.getvalue()


def guide_oligo_ordering_csv(result: DesignResult) -> str:
    """Return only the two cloning oligos derived from the selected guide."""
    if not result.sequence_complete:
        raise OrderingError("The design is not sequence-complete, so guide oligos are withheld.")
    return guide_oligos_csv(
        f"{_safe_stem(result)}_selected_guide",
        result.top_guide.spacer,
    )


def _zip_member(name: str, content: str) -> tuple[str, bytes]:
    return name, content.encode("utf-8")


def ordering_package_members(result: DesignResult) -> tuple[tuple[str, bytes], ...]:
    """Return the fixed, user-facing ordering package contents."""
    _require_ordering_ready(result)
    stem = _safe_stem(result)
    return (
        _zip_member(f"{stem}_twist_sequences.csv", twist_sequences_csv(result)),
        _zip_member(f"{stem}_guide_oligos.csv", guide_oligo_ordering_csv(result)),
        _zip_member(f"{stem}_genotyping_primers.csv", genotyping_primers_csv(result)),
        _zip_member(
            f"{stem}_assembled_plasmid.gb",
            assembled_plasmid_genbank(result, record_date=ORDERING_GENBANK_DATE),
        ),
        _zip_member(
            f"{stem}_wild_type_locus.gb",
            locus_context_genbank(result, "wild_type", record_date=ORDERING_GENBANK_DATE),
        ),
        _zip_member(
            f"{stem}_edited_locus.gb",
            locus_context_genbank(result, "edited", record_date=ORDERING_GENBANK_DATE),
        ),
    )


def ordering_package_zip(result: DesignResult) -> bytes:
    """Build a deterministic ZIP containing the six supported output files."""
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        for member_name, content in ordering_package_members(result):
            info = zipfile.ZipInfo(member_name, date_time=(1980, 1, 1, 0, 0, 0))
            info.compress_type = zipfile.ZIP_DEFLATED
            info.external_attr = 0o600 << 16
            archive.writestr(info, content)
    return output.getvalue()
