from __future__ import annotations

import csv
from datetime import datetime, timezone
import io
import json
from typing import Any

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqRecord import SeqRecord

from .models import DesignResult, GuideCandidate, HomologyArm
from .sequence import fasta_record, wrap_sequence


def guides_csv(result: DesignResult) -> str:
    output = io.StringIO()
    fields = [
        "rank",
        "spacer_5to3",
        "pam",
        "target_with_pam_5to3",
        "chromosome_strand",
        "target_interval_1based",
        "nick_boundary0",
        "distance_to_insertion_bp",
        "gc_percent",
        "poly_t",
        "activity_heuristic",
        "target_destroyed",
        "pam_destroyed",
        "initial_longest_retained_segment_nt",
        "final_longest_retained_segment_nt",
        "blocking_mutation_required_after_final_design",
        "blocking_mutation_note",
        "rationale",
    ]
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    for guide in result.guides:
        writer.writerow(
            {
                "rank": guide.rank,
                "spacer_5to3": guide.spacer,
                "pam": guide.pam,
                "target_with_pam_5to3": guide.target_with_pam,
                "chromosome_strand": guide.chromosome_strand,
                "target_interval_1based": guide.target_interval_1based,
                "nick_boundary0": guide.nick_boundary0,
                "distance_to_insertion_bp": guide.distance_to_insertion,
                "gc_percent": f"{guide.gc_percent:.1f}",
                "poly_t": guide.poly_t,
                "activity_heuristic": guide.activity_heuristic,
                "target_destroyed": guide.target_destroyed,
                "pam_destroyed": guide.pam_destroyed,
                "initial_longest_retained_segment_nt": guide.longest_retained_segment,
                "final_longest_retained_segment_nt": guide.final_longest_retained_segment,
                "blocking_mutation_required_after_final_design": guide.blocking_mutation_required,
                "blocking_mutation_note": guide.blocking_mutation_note,
                "rationale": guide.rationale,
            }
        )
    return output.getvalue()


def arms_fasta(result: DesignResult) -> str:
    records: list[str] = []
    for arm in (result.five_prime_arm, result.three_prime_arm):
        prefix = arm.name.replace("-", "_").replace(" ", "_")
        records.append(
            fasta_record(
                f"{result.gene_symbol}_{prefix}_gene_oriented_reference "
                f"{result.assembly} {arm.genomic_interval_1based}",
                arm.gene_oriented_sequence,
            )
        )
        if arm.final_gene_oriented_sequence != arm.gene_oriented_sequence:
            records.append(
                fasta_record(
                    f"{result.gene_symbol}_{prefix}_gene_oriented_FINAL_with_silent_mutations",
                    arm.final_gene_oriented_sequence,
                )
            )
        records.append(
            fasta_record(
                f"{result.gene_symbol}_{prefix}_chromosome_forward_reference "
                f"{result.assembly} {arm.genomic_interval_1based}",
                arm.chromosome_forward_sequence,
            )
        )

    payload = result.donor_payload.get("payload_sequence_5to3")
    if payload:
        records.append(
            fasta_record(
                f"Addgene_{result.backbone_addgene_id}_GGGGSAS_mNeonGreen_stop_payload",
                str(payload),
            )
        )
    uha = result.cloning_fragments.get("uha_synthesis_fragment_5to3")
    dha = result.cloning_fragments.get("dha_synthesis_fragment_5to3")
    donor = result.cloning_fragments.get("assembled_donor_insert_5to3")
    if uha:
        records.append(fasta_record(f"{result.gene_symbol}_UHA_synthesis_fragment_FINAL", str(uha)))
    if dha:
        records.append(fasta_record(f"{result.gene_symbol}_DHA_synthesis_fragment_FINAL", str(dha)))
    if donor:
        records.append(fasta_record(f"{result.gene_symbol}_assembled_donor_insert_simulation", str(donor)))
    plasmid = result.cloning_fragments.get("assembled_plasmid_5to3")
    if plasmid:
        records.append(
            fasta_record(
                f"{result.gene_symbol}_C_terminal_mNeonGreen_Addgene_169227_assembled_circular_plasmid",
                str(plasmid),
            )
        )
    if result.edited_cds_sequence:
        records.append(fasta_record(f"{result.gene_symbol}_predicted_fusion_CDS_without_terminal_stop", result.edited_cds_sequence))
    return "".join(records)


def assembled_plasmid_genbank(result: DesignResult) -> str:
    """Return an annotated GenBank representation of the simulated circular plasmid."""
    sequence = result.cloning_fragments.get("assembled_plasmid_5to3")
    features = result.cloning_fragments.get("assembled_plasmid_features")
    if not sequence or not features:
        raise ValueError("No complete assembled plasmid is available for GenBank export")

    record = SeqRecord(
        Seq(str(sequence)),
        id=f"{result.gene_symbol}_mNG_HDR",
        name=f"{result.gene_symbol}_mNG_HDR"[:16],
        description=(
            f"Simulated SapI Golden Gate assembly: {result.gene_symbol} C-terminal "
            f"GGGGSAS-mNeonGreen donor in Addgene #{result.backbone_addgene_id}"
        ),
    )
    record.annotations.update(
        {
            "molecule_type": "DNA",
            "topology": "circular",
            "data_file_division": "SYN",
            "date": datetime.now(timezone.utc).strftime("%d-%b-%Y").upper(),
            "source": "synthetic DNA construct",
            "organism": "synthetic DNA construct",
            "taxonomy": ["other sequences", "artificial sequences"],
            "comment": (
                "Computational prototype output. Generated from the uploaded Addgene #169227 "
                "SnapGene sequence and Bollen supplementary SapI adapter rules. Independently "
                "verify before experimental use."
            ),
        }
    )
    record.features.append(
        SeqFeature(
            FeatureLocation(0, len(record.seq), strand=1),
            type="source",
            qualifiers={
                "organism": ["synthetic DNA construct"],
                "mol_type": ["other DNA"],
                "label": ["assembled plasmid"],
            },
        )
    )
    for item in features:
        start0 = int(item["start0"])
        end0 = int(item["end0"])
        if not (0 <= start0 < end0 <= len(record.seq)):
            raise ValueError(
                f"Invalid assembled-plasmid feature coordinates: {item.get('label')} "
                f"{start0}:{end0}"
            )
        qualifiers: dict[str, list[str]] = {
            "label": [str(item.get("label", item.get("type", "feature")))],
        }
        note = item.get("note")
        if note:
            qualifiers["note"] = [str(note)]
        if item.get("type") == "CDS" and "mNeonGreen" in str(item.get("label", "")):
            qualifiers.update(
                {
                    "product": ["mNeonGreen"],
                    "codon_start": ["1"],
                    "transl_table": ["1"],
                    "translation": [str(result.donor_payload.get("tag_peptide", ""))],
                }
            )
        record.features.append(
            SeqFeature(
                FeatureLocation(start0, end0, strand=int(item.get("strand", 0)) or None),
                type=str(item.get("type", "misc_feature")),
                qualifiers=qualifiers,
            )
        )

    coordinate_map = result.cloning_fragments.get("assembly_coordinate_map", {})
    uha_start0 = coordinate_map.get("uha_start0")
    if isinstance(uha_start0, int):
        for mutation in result.five_prime_arm.mutations:
            position0 = uha_start0 + mutation.arm_position1 - 1
            record.features.append(
                SeqFeature(
                    FeatureLocation(position0, position0 + 1, strand=1),
                    type="variation",
                    qualifiers={
                        "label": [mutation.kind],
                        "replace": [mutation.alternate_base],
                        "note": [
                            f"{mutation.reference_base}>{mutation.alternate_base}; "
                            f"{mutation.original_codon}>{mutation.altered_codon}; "
                            f"{mutation.amino_acid}; {mutation.reason}"
                        ],
                    },
                )
            )

    record.features.sort(key=lambda feature: (int(feature.location.start), int(feature.location.end)))
    handle = io.StringIO()
    SeqIO.write(record, handle, "genbank")
    return handle.getvalue()


def design_json(result: DesignResult) -> str:
    return json.dumps(result.to_dict(), indent=2, sort_keys=True) + "\n"


def _indent(text: str, prefix: str = "    ") -> str:
    return "\n".join(prefix + line for line in text.splitlines())


def _boundary_description(boundary0: int) -> str:
    return f"between 1-based bases {boundary0:,} and {boundary0 + 1:,}"


def _guide_report(guide: GuideCandidate) -> list[str]:
    return [
        f"Rank {guide.rank}: 5'-{guide.spacer}-{guide.pam}-3'",
        f"  Target + PAM: {guide.target_with_pam}",
        f"  PAM-containing chromosome strand: {guide.chromosome_strand}",
        f"  Genomic target interval (1-based): {guide.target_interval_1based}",
        f"  Nominal nick boundary (0-based): {guide.nick_boundary0:,}",
        f"  Distance to insertion: {guide.distance_to_insertion} bp",
        f"  GC: {guide.gc_percent:.1f}% | TTTT: {'yes' if guide.poly_t else 'no'}",
        f"  Activity heuristic: {guide.activity_heuristic}",
        f"  Intended edit disrupts target: {'yes' if guide.target_destroyed else 'no'}",
        f"  Intended edit disrupts PAM: {'yes' if guide.pam_destroyed else 'no'}",
        f"  Longest original segment after the intended edit: {guide.longest_retained_segment} nt",
        f"  Longest original segment after any donor-protection edits: {guide.final_longest_retained_segment} nt",
        f"  Blocking mutation note: {guide.blocking_mutation_note or '(none)'}",
        f"  Rationale: {guide.rationale}",
    ]


def _arm_report(arm: HomologyArm) -> list[str]:
    lines = [
        arm.name,
        f"  Reference interval (1-based, chromosome-forward): {arm.genomic_interval_1based}",
        f"  Length: {arm.length} bp | GC: {arm.gc_percent:.1f}%",
        f"  Raw internal SapI sites: {len(arm.sapi_sites)}",
        f"  Final internal SapI sites: {len(arm.final_sapi_sites)}",
    ]
    if arm.mutations:
        lines.append("  Final silent mutations:")
        for mutation in arm.mutations:
            lines.append(
                f"    - {mutation.kind}: arm base {mutation.arm_position1}, "
                f"chr{arm.chromosome}:{mutation.genomic_position1:,} "
                f"gene-oriented {mutation.reference_base}>{mutation.alternate_base}; "
                f"{mutation.original_codon}>{mutation.altered_codon}; {mutation.amino_acid}."
            )
            lines.append(f"      Reason: {mutation.reason}")
    lines.extend(
        [
            "  Reference gene-oriented sequence (5'->3'):",
            _indent(wrap_sequence(arm.gene_oriented_sequence)),
        ]
    )
    if arm.final_gene_oriented_sequence != arm.gene_oriented_sequence:
        lines.extend(
            [
                "  FINAL gene-oriented arm sequence (5'->3'):",
                _indent(wrap_sequence(arm.final_gene_oriented_sequence)),
            ]
        )
    return lines


def design_report(result: DesignResult) -> str:
    top = result.top_guide if result.guides else None
    lines: list[str] = [
        "HDR TAG DESIGNER - COMPUTATIONAL DESIGN REPORT",
        "=" * 58,
        "",
        f"Status: {result.status}",
        f"Sequence-complete computational output: {'yes' if result.sequence_complete else 'no'}",
        f"Species / assembly: {result.species_label} / {result.assembly}",
        f"Gene: {result.gene_symbol} ({result.gene_id})",
        f"Transcript: {result.transcript_id}",
        f"Locus: chr{result.chromosome}, gene strand {result.gene_strand}",
        f"Tagging mode: {result.terminus}",
        f"Nuclease: {result.nuclease_mode}",
        f"Backbone: {result.backbone_name}, Addgene #{result.backbone_addgene_id}",
        f"Homology arms: {result.homology_arm_length} bp each",
        "Off-target analysis: not performed (requested)",
        f"Guide-ranking note: {result.guide_scoring_note}",
        "Custom donor backbone support: deferred until after this validation test",
        "",
        "EDIT DEFINITION",
        "-" * 58,
        f"Insertion boundary: chr{result.chromosome}:{_boundary_description(result.insertion_boundary0)}",
        (
            f"Removed genomic interval (1-based): chr{result.chromosome}:"
            f"{result.removed_genomic_interval_start0 + 1:,}-"
            f"{result.removed_genomic_interval_end0:,}"
            if result.removed_genomic_interval_end0 > result.removed_genomic_interval_start0
            else "Removed genomic interval: none"
        ),
        f"Removed gene-oriented sequence: {result.removed_sequence_gene_oriented or '(none)'}",
        f"Endogenous coding sequence before tag: {result.cds_length_without_stop} nt / {result.protein_length_aa} aa",
        "",
        "SELECTED GUIDE",
        "-" * 58,
    ]
    if top:
        lines.extend(_guide_report(top))
    else:
        lines.append("No candidate guide found.")

    lines.extend(["", "TOP GUIDE CANDIDATES", "-" * 58])
    for guide in result.guides[:10]:
        lines.extend(_guide_report(guide))
        lines.append("")

    lines.extend(["HOMOLOGY ARMS", "-" * 58])
    lines.extend(_arm_report(result.five_prime_arm))
    lines.append("")
    lines.extend(_arm_report(result.three_prime_arm))

    if result.donor_payload:
        lines.extend(["", "FIXED DONOR PAYLOAD", "-" * 58])
        lines.extend(
            [
                f"Name: {result.donor_payload.get('name')}",
                f"Payload length: {result.donor_payload.get('payload_length_nt')} nt",
                f"Linker: {result.donor_payload.get('linker_coding_sequence')} / {result.donor_payload.get('linker_peptide')}",
                f"mNeonGreen coding region: {result.donor_payload.get('tag_length_nt')} nt / {result.donor_payload.get('tag_length_aa')} aa",
                f"Payload stop codon: {result.donor_payload.get('stop_codon')}",
                "Payload sequence (5'->3'):",
                _indent(wrap_sequence(str(result.donor_payload.get('payload_sequence_5to3', '')))),
            ]
        )

    if result.cloning_fragments:
        lines.extend(["", "SAPI / GOLDEN GATE OUTPUT", "-" * 58])
        for key in (
            "uha_synthesis_fragment_5to3",
            "dha_synthesis_fragment_5to3",
            "assembled_donor_insert_5to3",
            "assembled_plasmid_5to3",
            "uha_synthesis_fragment_preview_5to3",
            "dha_synthesis_fragment_preview_5to3",
        ):
            sequence = result.cloning_fragments.get(key)
            if sequence:
                lines.append(f"{key} ({len(str(sequence))} nt):")
                lines.append(_indent(wrap_sequence(str(sequence))))
        overhangs = result.cloning_fragments.get("expected_sapi_overhangs")
        if overhangs:
            lines.append(f"Expected overhangs: {overhangs}")
        fragment_sites = result.cloning_fragments.get("synthesis_fragment_sapi_sites")
        if fragment_sites:
            lines.append("Verified SapI digest of synthesis fragments:")
            for fragment_name, sites in fragment_sites.items():
                observed = ", ".join(str(site.get("overhang_5to3")) for site in sites)
                lines.append(f"  - {fragment_name}: {observed}")
        backbone = result.cloning_fragments.get("uploaded_backbone")
        if backbone:
            lines.extend(
                [
                    "Uploaded fixed-backbone verification:",
                    f"  File: {backbone.get('snapgene_file')}",
                    f"  SHA-256: {backbone.get('snapgene_sha256')}",
                    f"  Original plasmid: {backbone.get('length_nt')} bp, {backbone.get('topology')}",
                    f"  SapI sites: {backbone.get('sapi_site_count')}",
                    f"  Backbone payload matches Bollen S2: {backbone.get('payload_matches_bollen_s2')}",
                    f"  Simulated final plasmid: {result.cloning_fragments.get('assembled_plasmid_length_nt')} bp, circular",
                    f"  SapI sites in final circular plasmid: {result.cloning_fragments.get('assembled_plasmid_sapi_site_count')}",
                ]
            )
            for site in backbone.get("sapi_sites", []):
                lines.append(
                    "  - "
                    f"{site.get('motif')} at {site.get('recognition_interval_1based')} "
                    f"({site.get('orientation')}), cut boundary {site.get('top_strand_cut_boundary0')}, "
                    f"overhang {site.get('overhang_5to3')}"
                )
        ligation_junctions = result.cloning_fragments.get("golden_gate_junctions")
        if ligation_junctions:
            lines.append("Verified ligation junctions:")
            for name, details in ligation_junctions.items():
                lines.append(
                    f"  - {name}: expected {details.get('overhang_5to3')}, "
                    f"observed {details.get('observed')}; window {details.get('window_5to3')}"
                )

    if result.primer_tail_templates:
        lines.extend(["", "PCR PRIMER 5' TAIL TEMPLATES", "-" * 58])
        for key, value in result.primer_tail_templates.items():
            lines.append(f"{key}: {value}")

    if result.fusion_protein_sequence:
        lines.extend(["", "PREDICTED FUSION", "-" * 58])
        lines.append(f"Fusion length: {result.fusion_protein_length_aa} aa")
        lines.append("Protein sequence:")
        lines.append(_indent("\n".join(
            result.fusion_protein_sequence[i:i + 60]
            for i in range(0, len(result.fusion_protein_sequence), 60)
        )))
        lines.append("Junction windows:")
        for key, value in result.junctions.items():
            lines.append(f"  {key}:")
            lines.append(_indent(wrap_sequence(value), "    "))

    lines.extend(["", "VALIDATION GATES", "-" * 58])
    for validation in result.validations:
        lines.append(f"[{validation['status']}] {validation['check']}: {validation['detail']}")

    lines.extend(["", "WARNINGS", "-" * 58])
    for warning in result.warnings:
        lines.append(f"- {warning}")

    lines.extend(["", "PROVENANCE", "-" * 58])
    for item in result.provenance:
        lines.append(f"- {item}")

    lines.extend(
        [
            "",
            "COORDINATE CONVENTIONS",
            "-" * 58,
            "Displayed genomic intervals are 1-based inclusive. JSON uses 0-based half-open intervals and 0-based boundaries.",
            "",
        ]
    )
    return "\n".join(lines)


def guide_rows(result: DesignResult) -> list[dict[str, Any]]:
    return [
        {
            "Rank": guide.rank,
            "Spacer (5'->3')": guide.spacer,
            "PAM": guide.pam,
            "Chr strand": guide.chromosome_strand,
            "Nick distance (bp)": guide.distance_to_insertion,
            "GC %": guide.gc_percent,
            "Activity heuristic": guide.activity_heuristic,
            "Edit disrupts target": guide.target_destroyed,
            "Retained after intended edit (nt)": guide.longest_retained_segment,
            "Retained after donor edits (nt)": guide.final_longest_retained_segment,
            "TTTT": guide.poly_t,
        }
        for guide in result.guides
    ]
