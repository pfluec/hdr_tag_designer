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
        "final_pam",
        "final_pam_destroyed",
        "final_target_after_point_mutations_5to3",
        "edited_target_region_display_5to3",
        "edited_target_region_full_5to3",
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
                "final_pam": guide.final_pam,
                "final_pam_destroyed": guide.final_pam_destroyed,
                "final_target_after_point_mutations_5to3": guide.final_target_with_pam_after_point_mutations,
                "edited_target_region_display_5to3": guide.edited_target_region_display,
                "edited_target_region_full_5to3": guide.edited_target_region_5to3,
                "initial_longest_retained_segment_nt": guide.longest_retained_segment,
                "final_longest_retained_segment_nt": guide.final_longest_retained_segment,
                "blocking_mutation_required_after_final_design": guide.blocking_mutation_required,
                "blocking_mutation_note": guide.blocking_mutation_note,
                "rationale": guide.rationale,
            }
        )
    return output.getvalue()


def genotyping_primers_csv(result: DesignResult) -> str:
    """Return one row per primer, with assay-level expected product metadata."""
    output = io.StringIO()
    fields = [
        "assay",
        "assay_status",
        "used_in_assays",
        "primer_role",
        "primer_name",
        "sequence_5to3",
        "orientation",
        "reference_sequence_strand",
        "source",
        "reusable_payload_primer",
        "outside_homology_arm",
        "genomic_interval_1based",
        "payload_interval_1based",
        "distance_from_5prime_junction_nt",
        "distance_from_3prime_junction_nt",
        "length_nt",
        "tm_c",
        "gc_percent",
        "hairpin_tm_c",
        "homodimer_tm_c",
        "pair_heterodimer_tm_c",
        "pair_tm_difference_c",
        "expected_product_size_bp",
        "expected_wild_type_product_size_bp",
        "expected_edited_product_size_bp",
        "expected_amplicon_sequence_5to3",
        "expected_wild_type_amplicon_sequence_5to3",
        "expected_edited_amplicon_sequence_5to3",
        "reason",
    ]
    writer = csv.DictWriter(output, fieldnames=fields, lineterminator="\n")
    writer.writeheader()
    assays = result.genotyping_primers.get("assays", {})
    unique_primers = (
        ("wild_type_locus", "forward_primer", "wild_type_locus;five_prime_junction"),
        ("wild_type_locus", "reverse_primer", "wild_type_locus;three_prime_junction"),
        ("five_prime_junction", "reverse_primer", "five_prime_junction"),
        ("three_prime_junction", "forward_primer", "three_prime_junction"),
    )
    for assay_name, role, used_in_assays in unique_primers:
        assay = assays.get(assay_name, {})
        primer = assay.get(role)
        if not primer:
            continue
        writer.writerow(
            {
                "assay": assay_name,
                "assay_status": assay.get("status", ""),
                "used_in_assays": used_in_assays,
                "primer_role": role.replace("_primer", ""),
                "primer_name": primer.get("name", ""),
                "sequence_5to3": primer.get("sequence_5to3", ""),
                "orientation": primer.get("orientation", ""),
                "reference_sequence_strand": primer.get("reference_sequence_strand", ""),
                "source": primer.get("source", ""),
                "reusable_payload_primer": primer.get("source") == "payload",
                "outside_homology_arm": primer.get("outside_homology_arm", False),
                "genomic_interval_1based": primer.get("genomic_interval_1based", ""),
                "payload_interval_1based": primer.get("payload_interval_1based", ""),
                "distance_from_5prime_junction_nt": primer.get("distance_from_5prime_junction_nt", ""),
                "distance_from_3prime_junction_nt": primer.get("distance_from_3prime_junction_nt", ""),
                "length_nt": primer.get("length_nt", ""),
                "tm_c": primer.get("tm_c", ""),
                "gc_percent": primer.get("gc_percent", ""),
                "hairpin_tm_c": primer.get("hairpin_tm_c", ""),
                "homodimer_tm_c": primer.get("homodimer_tm_c", ""),
                "pair_heterodimer_tm_c": assay.get("heterodimer_tm_c", ""),
                "pair_tm_difference_c": assay.get("tm_difference_c", ""),
                "expected_product_size_bp": assay.get("product_size_bp", ""),
                "expected_wild_type_product_size_bp": assay.get("expected_wild_type_product_size_bp", ""),
                "expected_edited_product_size_bp": assay.get("expected_edited_product_size_bp", ""),
                "expected_amplicon_sequence_5to3": assay.get("amplicon_sequence_5to3", ""),
                "expected_wild_type_amplicon_sequence_5to3": assay.get("expected_wild_type_amplicon_sequence_5to3", ""),
                "expected_edited_amplicon_sequence_5to3": assay.get("expected_edited_amplicon_sequence_5to3", ""),
                "reason": assay.get("reason", ""),
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
                f"Addgene_{result.backbone_addgene_id}_{result.terminus.replace('-', '_')}_payload",
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
                f"{result.gene_symbol}_{result.terminus.replace('-', '_')}_"
                f"{str(result.donor_payload.get('tag_name', 'tag')).replace(' ', '_')}_"
                f"Addgene_{result.backbone_addgene_id}_assembled_circular_plasmid",
                str(plasmid),
            )
        )
    if result.edited_cds_sequence:
        records.append(fasta_record(f"{result.gene_symbol}_predicted_fusion_CDS_without_terminal_stop", result.edited_cds_sequence))
    for assay_name, assay in result.genotyping_primers.get("assays", {}).items():
        amplicon = assay.get("amplicon_sequence_5to3")
        if amplicon:
            records.append(
                fasta_record(
                    f"{result.gene_symbol}_{assay_name}_expected_amplicon",
                    str(amplicon),
                )
            )
        if assay_name == "wild_type_locus":
            edited_amplicon = assay.get("expected_edited_amplicon_sequence_5to3")
            if edited_amplicon:
                records.append(
                    fasta_record(
                        f"{result.gene_symbol}_wild_type_locus_expected_edited_amplicon",
                        str(edited_amplicon),
                    )
                )
    for context_name, context in result.locus_contexts.items():
        if not isinstance(context, dict) or not context.get("sequence_5to3"):
            continue
        records.append(
            fasta_record(
                f"{result.gene_symbol}_{context_name}_locus_with_300bp_external_flanks_gene_oriented",
                str(context["sequence_5to3"]),
            )
        )
    return "".join(records)


def locus_context_genbank(
    result: DesignResult,
    context_name: str,
    *,
    record_date: str | None = None,
) -> str:
    """Return one annotated linear WT or edited locus context as GenBank."""
    context = result.locus_contexts.get(context_name)
    if not isinstance(context, dict) or not context.get("sequence_5to3"):
        raise ValueError(f"No {context_name} locus context is available")
    sequence = str(context["sequence_5to3"])
    suffix = "WT" if context_name == "wild_type" else "EDITED"
    record = SeqRecord(
        Seq(sequence),
        id=f"{result.gene_symbol}_{suffix}_locus",
        name=f"{result.gene_symbol}_{suffix}"[:16],
        description=(
            f"{result.gene_symbol} {suffix.lower()} gene-oriented locus context with "
            "300-bp sequence beyond each homology arm"
        ),
    )
    record.annotations.update(
        {
            "molecule_type": "DNA",
            "topology": "linear",
            "data_file_division": "SYN",
            "date": record_date or datetime.now(timezone.utc).strftime("%d-%b-%Y").upper(),
            "source": f"{result.species_label} reference and computational donor design",
            "organism": result.species_label,
            "comment": (
                f"Gene-oriented 5-prime to 3-prime context on {result.assembly}. "
                "Primer annotations are computational candidates; confirm genome-wide "
                "specificity before ordering."
            ),
        }
    )
    record.features.append(
        SeqFeature(
            FeatureLocation(0, len(sequence), strand=1),
            type="source",
            qualifiers={"label": [f"{suffix.lower()} locus context"]},
        )
    )
    for item in context.get("features", []):
        start0 = int(item["start0"])
        end0 = int(item["end0"])
        if not 0 <= start0 < end0 <= len(sequence):
            raise ValueError(
                f"Invalid {context_name} locus feature coordinates: "
                f"{item.get('label')} {start0}:{end0}"
            )
        qualifiers = {"label": [str(item.get("label", item.get("type", "feature")))]}
        if item.get("note"):
            qualifiers["note"] = [str(item["note"])]
        record.features.append(
            SeqFeature(
                FeatureLocation(start0, end0, strand=int(item.get("strand", 1))),
                type=str(item.get("type", "misc_feature")),
                qualifiers=qualifiers,
            )
        )
    record.features.sort(key=lambda item: (int(item.location.start), int(item.location.end)))
    handle = io.StringIO()
    SeqIO.write(record, handle, "genbank")
    return handle.getvalue()


def assembled_plasmid_genbank(
    result: DesignResult,
    *,
    record_date: str | None = None,
) -> str:
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
            f"Simulated SapI Golden Gate assembly: {result.gene_symbol} {result.terminus} "
            f"{result.donor_payload.get('tag_name', 'tag')} donor in "
            f"{result.backbone_name}"
        ),
    )
    record.annotations.update(
        {
            "molecule_type": "DNA",
            "topology": "circular",
            "data_file_division": "SYN",
            "date": record_date or datetime.now(timezone.utc).strftime("%d-%b-%Y").upper(),
            "source": "synthetic DNA construct",
            "organism": "synthetic DNA construct",
            "taxonomy": ["other sequences", "artificial sequences"],
            "comment": (
                f"Computational prototype output. Generated from the uploaded Addgene "
                f"#{result.backbone_addgene_id} SnapGene sequence and Bollen supplementary "
                "SapI adapter rules. Independently "
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
        tag_name = str(result.donor_payload.get("tag_name", "mNeonGreen"))
        if item.get("type") == "CDS" and str(item.get("label", "")) in {
            tag_name,
            f"{tag_name} + stop",
        }:
            qualifiers.update(
                {
                    "product": [tag_name],
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
    for arm, coordinate_key in (
        (result.five_prime_arm, "uha_start0"),
        (result.three_prime_arm, "dha_start0"),
    ):
        arm_start0 = coordinate_map.get(coordinate_key)
        if not isinstance(arm_start0, int):
            continue
        for mutation in arm.mutations:
            position0 = arm_start0 + mutation.arm_position1 - 1
            metadata = [
                f"{mutation.reference_base}>{mutation.alternate_base}",
                f"{mutation.original_codon}>{mutation.altered_codon}",
                mutation.protein_consequence or mutation.amino_acid,
            ]
            if mutation.pam_before or mutation.pam_after:
                metadata.append(f"PAM {mutation.pam_before}>{mutation.pam_after}")
            if mutation.longest_retained_before is not None:
                metadata.append(
                    "longest retained target segment "
                    f"{mutation.longest_retained_before}>{mutation.longest_retained_after} nt"
                )
            metadata.append(mutation.reason)
            record.features.append(
                SeqFeature(
                    FeatureLocation(position0, position0 + 1, strand=1),
                    type="variation",
                    qualifiers={
                        "label": [mutation.kind],
                        "replace": [mutation.alternate_base],
                        "note": ["; ".join(part for part in metadata if part)],
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


def _sapi_site_genomic_interval(
    result: DesignResult,
    arm: HomologyArm,
    position0: int,
    motif_length: int,
) -> str:
    if result.gene_strand == "+":
        genomic_positions0 = range(
            arm.genomic_start0 + position0,
            arm.genomic_start0 + position0 + motif_length,
        )
    else:
        first0 = arm.genomic_end0 - 1 - position0
        genomic_positions0 = range(first0 - motif_length + 1, first0 + 1)
    return (
        f"chr{arm.chromosome}:{min(genomic_positions0) + 1:,}-"
        f"{max(genomic_positions0) + 1:,}"
    )


def sapi_qc_rows(result: DesignResult) -> list[dict[str, Any]]:
    """Return one quality-control row for every SapI site originally found in an arm."""
    rows: list[dict[str, Any]] = []
    for arm in (result.five_prime_arm, result.three_prime_arm):
        final_site_keys = {
            (str(site["motif"]), int(site["position0"]))
            for site in arm.final_sapi_sites
        }
        raw_site_keys = {
            (str(site["motif"]), int(site["position0"]))
            for site in arm.sapi_sites
        }
        for site_number, site in enumerate(arm.sapi_sites, start=1):
            motif = str(site["motif"])
            position0 = int(site["position0"])
            direct_mutations = [
                mutation
                for mutation in arm.mutations
                if mutation.kind == "SapI domestication"
                and position0 <= mutation.arm_position1 - 1 < position0 + len(motif)
            ]
            related_keys = {
                (mutation.original_codon, mutation.altered_codon, mutation.reason)
                for mutation in direct_mutations
            }
            related_mutations = [
                mutation
                for mutation in arm.mutations
                if mutation.kind == "SapI domestication"
                and (
                    mutation in direct_mutations
                    or (
                        (mutation.original_codon, mutation.altered_codon, mutation.reason)
                        in related_keys
                        and any(
                            abs(mutation.arm_position1 - direct.arm_position1) <= 2
                            for direct in direct_mutations
                        )
                    )
                )
            ]
            resolved = (motif, position0) not in final_site_keys
            mutation_text = "; ".join(
                f"arm base {mutation.arm_position1}, "
                f"chr{arm.chromosome}:{mutation.genomic_position1:,} "
                f"{mutation.reference_base}>{mutation.alternate_base}"
                for mutation in related_mutations
            ) or "None"
            codon_changes = sorted(
                {
                    f"{mutation.original_codon}>{mutation.altered_codon}"
                    for mutation in related_mutations
                    if mutation.original_codon and mutation.altered_codon
                }
            )
            protein_consequences = sorted(
                {
                    mutation.protein_consequence or mutation.amino_acid
                    for mutation in related_mutations
                    if mutation.protein_consequence or mutation.amino_acid
                }
            )
            rows.append(
                {
                    "Arm": arm.name,
                    "Site": site_number,
                    "Motif": motif,
                    "Arm interval": f"{position0 + 1}-{position0 + len(motif)}",
                    "Genomic interval": _sapi_site_genomic_interval(
                        result, arm, position0, len(motif)
                    ),
                    "Status": "Resolved" if resolved else "Unresolved - design blocked",
                    "Mutation(s)": mutation_text,
                    "Codon change": ", ".join(codon_changes) or "Not applicable",
                    "Protein consequence": (
                        ", ".join(protein_consequences) or "Not established"
                    ),
                    "Selection reason": (
                        related_mutations[0].reason
                        if related_mutations
                        else "No verified synonymous correction was available."
                    ),
                }
            )
        for motif, position0 in sorted(final_site_keys - raw_site_keys, key=lambda item: item[1]):
            rows.append(
                {
                    "Arm": arm.name,
                    "Site": "new",
                    "Motif": motif,
                    "Arm interval": f"{position0 + 1}-{position0 + len(motif)}",
                    "Genomic interval": _sapi_site_genomic_interval(
                        result, arm, position0, len(motif)
                    ),
                    "Status": "New site introduced - design blocked",
                    "Mutation(s)": "See complete arm mutation list",
                    "Codon change": "Not released",
                    "Protein consequence": "Not released",
                    "Selection reason": "A SapI site is present only in the final arm sequence.",
                }
            )
    return rows


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
        f"  Final donor PAM: {guide.final_pam or guide.pam}",
        f"  Final donor disrupts PAM: {'yes' if guide.final_pam_destroyed else 'no'}",
        f"  Reference guide-binding region + PAM (5'->3'): {guide.target_with_pam}",
        f"  Final target after donor point mutations (5'->3'): {guide.final_target_with_pam_after_point_mutations or '(not reconstructed)'}",
        f"  Actual edited target region (compact, 5'->3'): {guide.edited_target_region_display or '(not reconstructed)'}",
        f"  Deleted target-region bases: {guide.edited_target_deleted_bases or '(none)'}",
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
        lines.append("  Final automatic sequence changes:")
        for mutation in arm.mutations:
            consequence = (
                f"{mutation.original_codon}>{mutation.altered_codon}; {mutation.amino_acid}"
                if mutation.original_codon and mutation.altered_codon
                else mutation.protein_consequence or "non-coding"
            )
            lines.append(
                f"    - {mutation.kind}: arm base {mutation.arm_position1}, "
                f"chr{arm.chromosome}:{mutation.genomic_position1:,} "
                f"gene-oriented {mutation.reference_base}>{mutation.alternate_base}; "
                f"{consequence}."
            )
            lines.append(f"      Reason: {mutation.reason}")
            if mutation.protein_consequence:
                lines.append(f"      Protein consequence: {mutation.protein_consequence}")
            if mutation.pam_before or mutation.pam_after:
                lines.append(
                    f"      PAM before/after: {mutation.pam_before or '(not applicable)'} -> "
                    f"{mutation.pam_after or '(not applicable)'}"
                )
            if mutation.longest_retained_before is not None:
                lines.append(
                    "      Longest retained target segment before/after: "
                    f"{mutation.longest_retained_before} -> "
                    f"{mutation.longest_retained_after} nt"
                )
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
    backbone_identifier = (
        f"Addgene #{result.backbone_addgene_id}"
        if result.backbone_addgene_id != "custom"
        else "user-supplied backbone"
    )
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
        f"Backbone: {result.backbone_name}, {backbone_identifier}",
        f"Homology arms: {result.homology_arm_length} bp each",
        "Off-target analysis: not performed (requested)",
        f"Guide-ranking note: {result.guide_scoring_note}",
        f"Custom donor backbone support: {'enabled' if result.custom_backbones_supported else 'disabled'}",
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
        f"Native coding sequence without terminal stop: {result.cds_length_without_stop} nt / {result.protein_length_aa} aa",
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

    sapi_rows = sapi_qc_rows(result)
    sapi_found = sum(
        len(arm.sapi_sites)
        for arm in (result.five_prime_arm, result.three_prime_arm)
    )
    sapi_remaining = sum(
        len(arm.final_sapi_sites)
        for arm in (result.five_prime_arm, result.three_prime_arm)
    )
    lines.extend(["SAPI ARM QUALITY CONTROL", "-" * 58])
    lines.append(f"Original SapI sites found in both arms: {sapi_found}")
    lines.append(f"Original SapI sites resolved: {sum(row['Status'] == 'Resolved' for row in sapi_rows)}")
    lines.append(f"SapI sites remaining in final arms: {sapi_remaining}")
    if not sapi_rows:
        lines.append("No GCTCTTC/GAAGAGC recognition motif was found in either arm.")
    for row in sapi_rows:
        lines.extend(
            [
                f"- {row['Arm']} site {row['Site']}: {row['Motif']} at arm bases "
                f"{row['Arm interval']} ({row['Genomic interval']})",
                f"  Status: {row['Status']}",
                f"  Mutation(s): {row['Mutation(s)']}",
                f"  Codon / protein: {row['Codon change']}; {row['Protein consequence']}",
                f"  Reason: {row['Selection reason']}",
            ]
        )

    lines.extend(["", "HOMOLOGY ARMS", "-" * 58])
    lines.extend(_arm_report(result.five_prime_arm))
    lines.append("")
    lines.extend(_arm_report(result.three_prime_arm))

    if result.donor_payload:
        lines.extend(["", "DONOR PAYLOAD", "-" * 58])
        lines.extend(
            [
                f"Name: {result.donor_payload.get('name')}",
                f"Interpretation: {result.donor_payload.get('payload_kind')}",
                f"Interpretation warning: {result.donor_payload.get('payload_warning') or '(none)'}",
                f"Payload length: {result.donor_payload.get('payload_length_nt')} nt",
                f"Linker: {result.donor_payload.get('linker_coding_sequence')} / {result.donor_payload.get('linker_peptide')}",
                f"{result.donor_payload.get('tag_name', 'Tag')} coding region: "
                f"{result.donor_payload.get('tag_length_nt')} nt / "
                f"{result.donor_payload.get('tag_length_aa')} aa",
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
                    "Uploaded backbone verification:",
                    f"  File: {backbone.get('snapgene_file')}",
                    f"  SHA-256: {backbone.get('snapgene_sha256')}",
                    f"  Original plasmid: {backbone.get('length_nt')} bp, {backbone.get('topology')}",
                    f"  SapI sites: {backbone.get('sapi_site_count')}",
                    f"  Backbone payload passed sequence checks: {backbone.get('payload_sequence_verified')}",
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

    if result.cloning_primers:
        lines.extend(["", "HOMOLOGY-ARM CLONING PRIMERS (SAPI / GOLDEN GATE)", "-" * 58])
        lines.append(f"Status: {result.cloning_primers.get('status')}")
        lines.append(f"Rule set: {result.cloning_primers.get('ruleset')}")
        for primer in result.cloning_primers.get("primers", {}).values():
            lines.extend(
                [
                    f"{primer.get('name')}:",
                    f"  5' tail: {primer.get('tail_sequence_5to3')}",
                    f"  Genomic annealing region: {primer.get('annealing_sequence_5to3')}",
                    f"  Complete primer (5'->3'): {primer.get('full_sequence_5to3')}",
                    f"  Annealing length/Tm/GC: {primer.get('annealing_length_nt')} nt / "
                    f"{primer.get('annealing_tm_c')} C / {primer.get('annealing_gc_percent')}%",
                    f"  Arm binding interval: {primer.get('arm_binding_interval_1based')}",
                ]
            )
        for warning in result.cloning_primers.get("warnings", []):
            lines.append(f"Warning: {warning}")

    if result.genotyping_primers:
        lines.extend(["", "GENOTYPING PCR ASSAYS", "-" * 58])
        lines.append(f"Overall status: {result.genotyping_primers.get('status')}")
        lines.append(f"Rule set: {result.genotyping_primers.get('ruleset')}")
        lines.append(
            "Reference: "
            f"{result.genotyping_primers.get('assembly')} chr"
            f"{result.genotyping_primers.get('chromosome')}"
        )
        lines.append(
            "Minimum payload-primer distance from the tested junction: "
            f"{result.genotyping_primers.get('payload_junction_standoff_nt')} bp"
        )
        for assay_name, assay in result.genotyping_primers.get("assays", {}).items():
            lines.append(f"{assay_name}: {assay.get('status')}")
            if assay.get("reason"):
                lines.append(f"  Reason: {assay.get('reason')}")
                continue
            lines.append(f"  Expected product: {assay.get('product_size_bp')} bp")
            if assay.get("expected_wild_type_product_size_bp"):
                lines.append(
                    "  Expected WT / edited products: "
                    f"{assay.get('expected_wild_type_product_size_bp')} / "
                    f"{assay.get('expected_edited_product_size_bp')} bp"
                )
            for role in ("forward_primer", "reverse_primer"):
                primer = assay.get(role, {})
                if not primer:
                    continue
                location = primer.get("genomic_interval_1based") or primer.get("payload_interval_1based")
                lines.append(
                    f"  {role}: 5'-{primer.get('sequence_5to3')}-3'; "
                    f"{primer.get('length_nt')} nt; Tm {primer.get('tm_c')} C; "
                    f"GC {primer.get('gc_percent')}%; {primer.get('source')}; "
                    f"location {location}"
                )
                if primer.get("source") == "payload":
                    lines.append(
                        "    Payload distances from 5'/3' junctions: "
                        f"{primer.get('distance_from_5prime_junction_nt')} / "
                        f"{primer.get('distance_from_3prime_junction_nt')} bp"
                    )
            for allele, behavior in assay.get("allele_behavior", {}).items():
                lines.append(f"  {allele}: {behavior}")
            if assay.get("expected_wild_type_amplicon_sequence_5to3"):
                lines.append("  Expected WT amplicon sequence (5'->3'):")
                lines.append(
                    _indent(
                        wrap_sequence(
                            str(assay["expected_wild_type_amplicon_sequence_5to3"])
                        )
                    )
                )
                lines.append("  Expected edited amplicon sequence (5'->3'):")
                lines.append(
                    _indent(
                        wrap_sequence(
                            str(assay["expected_edited_amplicon_sequence_5to3"])
                        )
                    )
                )
            else:
                lines.append("  Expected amplicon sequence (5'->3'):")
                lines.append(
                    _indent(wrap_sequence(str(assay.get("amplicon_sequence_5to3", ""))))
                )
        for warning in result.genotyping_primers.get("warnings", []):
            lines.append(f"Warning: {warning}")

    if result.locus_contexts:
        lines.extend(["", "WT AND EDITED LOCUS CONTEXTS", "-" * 58])
        lines.append(
            "Orientation: gene-oriented 5'->3'; both records extend 300 bp beyond each homology arm."
        )
        for context_name in ("wild_type", "edited"):
            context = result.locus_contexts.get(context_name, {})
            lines.append(f"{context_name} locus ({context.get('length_nt')} bp):")
            lines.append("  Annotations:")
            for item in context.get("features", []):
                strand = "+" if item.get("strand", 1) == 1 else "-"
                lines.append(
                    f"    - {item.get('label')}: {int(item.get('start0', 0)) + 1}-"
                    f"{item.get('end0')} ({strand}); {item.get('type')}"
                )
            lines.append("  Sequence (5'->3'):")
            lines.append(_indent(wrap_sequence(str(context.get("sequence_5to3", "")))))

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
