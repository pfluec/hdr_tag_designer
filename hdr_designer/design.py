from __future__ import annotations

from dataclasses import replace

from .backbones import (
    BACKBONE_ADDGENE_ID,
    BACKBONE_NAME,
    CTERM_DHA_FORWARD_PRIMER_TAIL,
    CTERM_DHA_REVERSE_PRIMER_TAIL_PREFIX,
    CTERM_UHA_FORWARD_PRIMER_TAIL_PREFIX,
    CTERM_UHA_REVERSE_PRIMER_TAIL,
    NTERM_DHA_PREFIX,
    NTERM_DHA_SUFFIX,
    NTERM_UHA_PREFIX,
    NTERM_UHA_SUFFIX,
    SAPI_OVERHANGS,
    SAPI_RECOGNITION_MOTIFS,
    c_terminal_synthesis_fragments,
    payload_metadata,
)
from .ensembl import EnsemblClient, SPECIES
from .fixtures import (
    TUBB5_CDS_START_CDNA0,
    TUBB5_CHROMOSOME,
    TUBB5_CURRENT_CANONICAL_VERSION_AT_BUILD,
    TUBB5_FIXTURE_VERSION,
    TUBB5_GENE_ID,
    TUBB5_INSERTION_BOUNDARY0,
    TUBB5_INSERTION_CDNA0,
    TUBB5_PROTEIN_LENGTH_AA,
    TUBB5_SPECIES,
    TUBB5_STOP_CODON,
    TUBB5_STOP_END0,
    TUBB5_STOP_START0,
    TUBB5_STRAND,
    TUBB5_TRANSCRIPT_ID,
    load_tubb5_cdna,
)
from .guides import GUIDE_SAFETY_CUTOFF_NT, enumerate_spcas9_guides
from .models import (
    DesignResult,
    Exon,
    GuideCandidate,
    HomologyArm,
    SequenceMutation,
    TranscriptRecord,
)
from .sequence import CODON_TABLE, find_motif_positions, gc_percent, reverse_complement, translate


class DesignError(RuntimeError):
    pass


def _ordered_exons(record: TranscriptRecord) -> list[Exon]:
    return sorted(record.exons, key=lambda exon: exon.start1, reverse=record.strand == -1)


def build_transcript_genome_map(record: TranscriptRecord) -> list[int]:
    """Map each cDNA base index to a 0-based genomic base coordinate."""
    mapping: list[int] = []
    for exon in _ordered_exons(record):
        if record.strand == 1:
            mapping.extend(range(exon.start1 - 1, exon.end1))
        else:
            mapping.extend(range(exon.end1 - 1, exon.start1 - 2, -1))
    if len(mapping) != len(record.cdna):
        raise DesignError(
            "Transcript/exon length mismatch: Ensembl returned a cDNA of "
            f"{len(record.cdna)} nt but exon coordinates total {len(mapping)} nt."
        )
    return mapping


def transcript_boundary_to_genome0(mapping: list[int], boundary_index: int, strand: int) -> int:
    if not 0 <= boundary_index <= len(mapping):
        raise DesignError(f"Transcript boundary {boundary_index} is outside the cDNA")
    if boundary_index == 0:
        return mapping[0] if strand == 1 else mapping[0] + 1
    if boundary_index == len(mapping):
        return mapping[-1] + 1 if strand == 1 else mapping[-1]
    left = mapping[boundary_index - 1]
    right = mapping[boundary_index]
    if abs(left - right) != 1:
        raise DesignError(
            "The requested insertion boundary lies across an exon junction. "
            "This prototype releases only insertion sites within one genomic exon."
        )
    return max(left, right)


def _find_unique_cds(cdna: str, cds: str) -> int:
    first = cdna.find(cds)
    if first < 0:
        raise DesignError("The Ensembl CDS could not be located within the transcript cDNA")
    if cdna.find(cds, first + 1) >= 0:
        raise DesignError("The Ensembl CDS occurs more than once in the transcript cDNA")
    return first


def _arm_site_records(sequence: str) -> list[dict[str, int | str]]:
    sites: list[dict[str, int | str]] = []
    for motif, positions in find_motif_positions(sequence, SAPI_RECOGNITION_MOTIFS).items():
        for position in positions:
            sites.append({"motif": motif, "position0": position, "position1": position + 1})
    return sorted(sites, key=lambda item: (int(item["position0"]), str(item["motif"])))


def _make_arm(
    *,
    name: str,
    chromosome: str,
    start0: int,
    end0: int,
    chromosome_forward_sequence: str,
    gene_strand: int,
) -> HomologyArm:
    gene_oriented = (
        chromosome_forward_sequence
        if gene_strand == 1
        else reverse_complement(chromosome_forward_sequence)
    )
    sites = _arm_site_records(gene_oriented)
    return HomologyArm(
        name=name,
        length=len(gene_oriented),
        chromosome=chromosome,
        genomic_start0=start0,
        genomic_end0=end0,
        gene_oriented_sequence=gene_oriented,
        chromosome_forward_sequence=chromosome_forward_sequence,
        gc_percent=gc_percent(gene_oriented),
        sapi_sites=sites,
        final_sapi_sites=sites,
    )


def _arm_index_to_genomic_position1(arm: HomologyArm, index0: int, strand: int) -> int:
    genomic0 = (
        arm.genomic_start0 + index0
        if strand == 1
        else arm.genomic_end0 - 1 - index0
    )
    return genomic0 + 1


def _apply_point_mutation(
    arm: HomologyArm,
    *,
    index0: int,
    alternate_base: str,
    strand: int,
    kind: str,
    transcript_position1: int | None,
    original_codon: str,
    altered_codon: str,
    amino_acid: str,
    reason: str,
) -> HomologyArm:
    sequence = list(arm.final_gene_oriented_sequence)
    if not 0 <= index0 < len(sequence):
        raise DesignError(f"Mutation index {index0} is outside {arm.name}")
    reference_base = sequence[index0]
    alternate_base = alternate_base.upper()
    if reference_base == alternate_base:
        raise DesignError("Proposed mutation does not change the sequence")
    sequence[index0] = alternate_base
    corrected = "".join(sequence)
    mutation = SequenceMutation(
        kind=kind,
        arm_name=arm.name,
        arm_position1=index0 + 1,
        genomic_position1=_arm_index_to_genomic_position1(arm, index0, strand),
        transcript_position1=transcript_position1,
        reference_base=reference_base,
        alternate_base=alternate_base,
        original_codon=original_codon,
        altered_codon=altered_codon,
        amino_acid=amino_acid,
        reason=reason,
    )
    note = (
        f"{kind}: gene-oriented {reference_base}>{alternate_base} at arm base {index0 + 1} "
        f"(chr{arm.chromosome}:{mutation.genomic_position1:,}); "
        f"{original_codon}->{altered_codon}, {amino_acid} unchanged."
    )
    return replace(
        arm,
        corrected_gene_oriented_sequence=corrected,
        final_sapi_sites=_arm_site_records(corrected),
        mutations=[*arm.mutations, mutation],
        correction_note=" ".join(part for part in [arm.correction_note, note] if part),
    )


def _tubb5_corrected_five_prime_arm(
    arm: HomologyArm, *, cdna: str, arm_cdna_start0: int
) -> tuple[HomologyArm, list[tuple[int, str]]]:
    """Domesticate the single internal SapI site in the Tubb5 5-prime arm."""
    changes: list[tuple[int, str]] = []

    # Internal SapI site at arm bases 396-402. The first base of the motif is
    # the wobble base of GAG (Glu 376); GAG -> GAA removes the site.
    sapi_index0 = 395
    if arm.gene_oriented_sequence[sapi_index0:sapi_index0 + 7] != "GCTCTTC":
        raise DesignError("Expected Tubb5 SapI site was not found at arm base 396")
    codon_start_arm0 = 393
    original = arm.gene_oriented_sequence[codon_start_arm0:codon_start_arm0 + 3]
    if original != "GAG":
        raise DesignError("Expected Tubb5 SapI-overlapping GAG codon was not found")
    arm = _apply_point_mutation(
        arm,
        index0=sapi_index0,
        alternate_base="A",
        strand=TUBB5_STRAND,
        kind="SapI domestication",
        transcript_position1=arm_cdna_start0 + sapi_index0 + 1,
        original_codon="GAG",
        altered_codon="GAA",
        amino_acid="Glu (E)",
        reason="Remove the internal GCTCTTC recognition site without changing Tubb5.",
    )
    changes.append((arm_cdna_start0 + sapi_index0, "A"))

    if arm.final_sapi_sites:
        raise DesignError("A SapI site remains in the corrected Tubb5 5-prime arm")
    return arm, changes


def _apply_cdna_changes_to_cds(
    cds_without_stop: str,
    *,
    cds_start_cdna0: int,
    changes: list[tuple[int, str]],
) -> str:
    edited = list(cds_without_stop)
    for cdna_position0, alternate in changes:
        cds_index0 = cdna_position0 - cds_start_cdna0
        if not 0 <= cds_index0 < len(edited):
            raise DesignError("Arm mutation could not be mapped into the CDS")
        edited[cds_index0] = alternate
    return "".join(edited)


def _preview_n_terminal_fragments(target: str, uha: str, dha: str) -> dict[str, str | int]:
    uha_fragment = NTERM_UHA_PREFIX + target + uha + NTERM_UHA_SUFFIX
    dha_fragment = NTERM_DHA_PREFIX + dha + target + NTERM_DHA_SUFFIX
    return {
        "uha_synthesis_fragment_preview_5to3": uha_fragment,
        "uha_synthesis_fragment_length_nt": len(uha_fragment),
        "dha_synthesis_fragment_preview_5to3": dha_fragment,
        "dha_synthesis_fragment_length_nt": len(dha_fragment),
        "note": "Locus-only preview; Addgene #169227 is a C-terminal backbone.",
    }


def _finalize_result(
    *,
    record: TranscriptRecord,
    terminus: str,
    arm_length: int,
    guide_window: int,
    insertion_boundary0: int,
    removed_start0: int,
    removed_end0: int,
    removed_gene_sequence: str,
    cds_without_stop: str,
    edited_cds_without_stop: str,
    five_prime_arm: HomologyArm,
    three_prime_arm: HomologyArm,
    guides: list[GuideCandidate],
    provenance: list[str],
    extra_warnings: list[str] | None = None,
) -> DesignResult:
    is_c_terminal = terminus.upper().startswith("C")
    top = guides[0] if guides else None
    final_sapi_count = len(five_prime_arm.final_sapi_sites) + len(three_prime_arm.final_sapi_sites)

    donor_payload: dict[str, object] = payload_metadata() if is_c_terminal else {}
    cloning_fragments: dict[str, object] = {}
    primer_tails: dict[str, str] = {}
    edited_cds = ""
    fusion_protein = ""
    junctions: dict[str, str] = {}

    guide_safe = bool(top and not top.blocking_mutation_required)
    can_release_cterm = bool(
        is_c_terminal and top and guide_safe and final_sapi_count == 0
    )
    if can_release_cterm and top:
        cloning_fragments = c_terminal_synthesis_fragments(
            target_with_pam=top.target_with_pam,
            uha=five_prime_arm.final_gene_oriented_sequence,
            dha=three_prime_arm.final_gene_oriented_sequence,
        )
        primer_tails = {
            "UHA_forward_5prime_tail": (
                CTERM_UHA_FORWARD_PRIMER_TAIL_PREFIX + top.target_with_pam
            ),
            "UHA_reverse_5prime_tail": CTERM_UHA_REVERSE_PRIMER_TAIL,
            "DHA_forward_5prime_tail": CTERM_DHA_FORWARD_PRIMER_TAIL,
            "DHA_reverse_5prime_tail": (
                CTERM_DHA_REVERSE_PRIMER_TAIL_PREFIX
                + reverse_complement(top.target_with_pam)
            ),
            "note": (
                "Append each tail to a separately designed locus-specific annealing sequence; "
                "this prototype does not design PCR annealing regions."
            ),
        }
        payload_coding = str(donor_payload["payload_coding_sequence"])
        edited_cds = edited_cds_without_stop + payload_coding
        fusion_protein = translate(edited_cds)
        payload_sequence = str(donor_payload["payload_sequence_5to3"])
        junctions = {
            "five_prime_junction_5to3": (
                five_prime_arm.final_gene_oriented_sequence[-60:] + payload_sequence[:60]
            ),
            "three_prime_junction_5to3": (
                payload_sequence[-60:] + three_prime_arm.final_gene_oriented_sequence[:60]
            ),
            "edited_locus_window_5to3": (
                five_prime_arm.final_gene_oriented_sequence[-100:]
                + payload_sequence
                + three_prime_arm.final_gene_oriented_sequence[:100]
            ),
        }
    elif not is_c_terminal and top:
        cloning_fragments = _preview_n_terminal_fragments(
            top.target_with_pam,
            five_prime_arm.final_gene_oriented_sequence,
            three_prime_arm.final_gene_oriented_sequence,
        )

    backbone_info = cloning_fragments.get("uploaded_backbone", {})
    backbone_verified = bool(
        is_c_terminal
        and backbone_info
        and backbone_info.get("addgene_id") == BACKBONE_ADDGENE_ID
        and backbone_info.get("length_nt") == 2768
        and backbone_info.get("topology") == "circular"
        and backbone_info.get("sapi_site_count") == 4
        and backbone_info.get("payload_matches_bollen_s2") is True
    )
    golden_gate_junctions = cloning_fragments.get("golden_gate_junctions", {})
    junctions_verified = bool(
        golden_gate_junctions
        and all(
            details.get("observed") == details.get("overhang_5to3")
            for details in golden_gate_junctions.values()
        )
    )
    plasmid_assembly_verified = bool(
        is_c_terminal
        and cloning_fragments.get("assembled_plasmid_length_nt")
        and cloning_fragments.get("assembled_plasmid_topology") == "circular"
        and cloning_fragments.get("assembled_plasmid_sapi_site_count") == 0
        and junctions_verified
    )

    validations: list[dict[str, str]] = [
        {
            "check": "Reference coding frame",
            "status": "PASS" if len(cds_without_stop) % 3 == 0 else "FAIL",
            "detail": f"Coding sequence before the tag is {len(cds_without_stop)} nt ({len(cds_without_stop) // 3} aa).",
        },
        {
            "check": "SpCas9-NGG candidates",
            "status": "PASS" if guides else "FAIL",
            "detail": (
                f"Found {len(guides)} candidate(s) with a nominal nick within +/-{guide_window} bp."
                if guides
                else "No candidate was found in the selected window."
            ),
        },
    ]
    if top:
        validations.extend(
            [
                {
                    "check": "Guide proximity",
                    "status": "PASS",
                    "detail": f"Selected guide nick is {top.distance_to_insertion} bp from the insertion boundary.",
                },
                {
                    "check": "Target disrupted by edit",
                    "status": "PASS" if top.target_destroyed else "REVIEW",
                    "detail": top.rationale,
                },
                {
                    "check": "Donor retargeting safeguard",
                    "status": "PASS" if guide_safe else "BLOCKED",
                    "detail": (
                        (
                            "The intended edit removes part of the NGG PAM; no additional "
                            "guide-blocking mutation is required. "
                            f"Longest surviving target segment: {top.final_longest_retained_segment} nt."
                        )
                        if top.pam_destroyed
                        else (
                            f"Longest uninterrupted original target segment after final donor edits: "
                            f"{top.final_longest_retained_segment} nt; protocol cutoff <= {GUIDE_SAFETY_CUTOFF_NT} nt."
                        )
                    ),
                },
            ]
        )
    validations.extend(
        [
            {
                "check": "Internal SapI sites after corrections",
                "status": "PASS" if final_sapi_count == 0 else "BLOCKED",
                "detail": (
                    "No GCTCTTC/GAAGAGC motifs remain in either final arm."
                    if final_sapi_count == 0
                    else f"{final_sapi_count} SapI recognition site(s) remain in the final arms."
                ),
            },
            {
                "check": "Endogenous C-terminal stop omitted from DHA",
                "status": (
                    "PASS"
                    if is_c_terminal and removed_gene_sequence in {"TAA", "TAG", "TGA"}
                    else ("N/A" if not is_c_terminal else "FAIL")
                ),
                "detail": (
                    f"Removed endogenous stop codon {removed_gene_sequence}; the fixed payload supplies {donor_payload.get('stop_codon', '')}."
                    if is_c_terminal
                    else "Not applicable to N-terminal tagging."
                ),
            },
            {
                "check": "Bollen SapI fragment architecture",
                "status": "PASS" if can_release_cterm else ("N/A" if not is_c_terminal else "BLOCKED"),
                "detail": (
                    "Exact S1/S3 C-terminal adapters and TAC/GGC/TGA/AAT overhangs were applied."
                    if can_release_cterm
                    else (
                        "N-terminal preview only; the selected fixed backbone is C-terminal."
                        if not is_c_terminal
                        else "Fragments are withheld until guide-retargeting and internal-SapI gates pass."
                    )
                ),
            },
            {
                "check": "Uploaded Addgene #169227 backbone",
                "status": "PASS" if backbone_verified else ("N/A" if not is_c_terminal else "BLOCKED"),
                "detail": (
                    f"Parsed {backbone_info.get('snapgene_file')} ({backbone_info.get('length_nt')} bp, "
                    f"{backbone_info.get('topology')}); four SapI sites yield TAC/GGC/TGA/AAT, "
                    "and the extracted 729-bp payload matches Bollen supplementary S2."
                    if backbone_verified
                    else (
                        "Not applicable to an N-terminal preview."
                        if not is_c_terminal
                        else "The uploaded fixed-backbone sequence did not pass structural verification."
                    )
                ),
            },
            {
                "check": "Full circular plasmid assembly",
                "status": "PASS" if plasmid_assembly_verified else ("N/A" if not is_c_terminal else "BLOCKED"),
                "detail": (
                    f"Reconstructed a {cloning_fragments.get('assembled_plasmid_length_nt')}-bp circular plasmid; "
                    "all four ligation junctions match and no SapI recognition site remains."
                    if plasmid_assembly_verified
                    else (
                        "Not applicable to an N-terminal preview."
                        if not is_c_terminal
                        else "No fully verified circular Golden Gate product was released."
                    )
                ),
            },
            {
                "check": "Fusion translation",
                "status": "PASS" if fusion_protein else ("N/A" if not is_c_terminal else "BLOCKED"),
                "detail": (
                    f"Predicted fusion is {len(fusion_protein)} aa; linker {donor_payload.get('linker_peptide', '')}; "
                    f"mNeonGreen {donor_payload.get('tag_length_aa', '')} aa; linker-plus-tag payload "
                    f"{donor_payload.get('payload_peptide_length_aa', '')} aa."
                    if fusion_protein
                    else "No final fusion translation released."
                ),
            },
        ]
    )

    sequence_complete = bool(
        can_release_cterm and fusion_protein and backbone_verified and plasmid_assembly_verified
    )
    status = "SEQUENCE-COMPLETE COMPUTATIONAL DESIGN" if sequence_complete else "REVIEW REQUIRED"
    warnings = list(extra_warnings or [])
    warnings.extend(
        [
            "No genome-wide off-target analysis is performed, as requested.",
            "Reference-genome sequence only; strain, cell-line, and clone-specific variants are not assessed.",
            "No experimental activity score is calculated. After the two primary Bollen priorities, ranking uses only GC/poly-T heuristics.",
            "Sequence-complete means the internal computational gates and uploaded-backbone assembly simulation passed; independently verify every sequence and plasmid junction before experimental use.",
            "Custom donor backbones and tags are deliberately not implemented in this version.",
        ]
    )
    if not is_c_terminal:
        warnings.append(
            "Addgene #169227 is C-terminal; N-terminal output is a locus-design preview, not a complete donor design."
        )

    return DesignResult(
        status=status,
        sequence_complete=sequence_complete,
        species_label=record.species.label,
        assembly=record.species.assembly,
        gene_symbol=record.gene_symbol,
        gene_id=record.gene_id,
        transcript_id=record.display_transcript_id,
        chromosome=record.chromosome,
        gene_strand="+" if record.strand == 1 else "-",
        terminus=terminus,
        nuclease_mode="ITPN with SpCas9 D10A nickase",
        backbone_name=BACKBONE_NAME,
        backbone_addgene_id=BACKBONE_ADDGENE_ID,
        insertion_boundary0=insertion_boundary0,
        removed_genomic_interval_start0=removed_start0,
        removed_genomic_interval_end0=removed_end0,
        removed_sequence_gene_oriented=removed_gene_sequence,
        cds_length_without_stop=len(cds_without_stop),
        protein_length_aa=len(cds_without_stop) // 3,
        guide_search_radius=guide_window,
        homology_arm_length=arm_length,
        guide_safety_cutoff_nt=GUIDE_SAFETY_CUTOFF_NT,
        guide_scoring_note=(
            "Rank: nick distance, target disruption, then basic spacer heuristics. "
            "No quantitative on-target or off-target model is run."
        ),
        guides=guides,
        five_prime_arm=five_prime_arm,
        three_prime_arm=three_prime_arm,
        donor_payload=donor_payload,
        cloning_fragments=cloning_fragments,
        primer_tail_templates=primer_tails,
        edited_cds_sequence=edited_cds,
        fusion_protein_sequence=fusion_protein,
        fusion_protein_length_aa=len(fusion_protein),
        junctions=junctions,
        validations=validations,
        warnings=warnings,
        provenance=provenance,
        custom_backbones_supported=False,
    )


def design_online(
    *,
    species_key: str,
    gene: str,
    transcript_id: str | None = None,
    terminus: str = "C-terminal",
    arm_length: int = 600,
    guide_window: int = 50,
    client: EnsemblClient | None = None,
) -> DesignResult:
    if species_key not in SPECIES:
        raise DesignError(f"Unsupported species: {species_key}")
    if arm_length < 100:
        raise DesignError("Homology arms shorter than 100 bp are not supported")
    client = client or EnsemblClient()
    record = client.transcript_record(SPECIES[species_key], gene, transcript_id)
    mapping = build_transcript_genome_map(record)
    cds_start = _find_unique_cds(record.cdna, record.cds)
    terminus_upper = terminus.upper()

    if terminus_upper.startswith("C"):
        stop = record.cds[-3:]
        if stop not in {"TAA", "TAG", "TGA"}:
            raise DesignError(f"Transcript CDS does not end in a standard stop codon: {stop}")
        insertion_t0 = cds_start + len(record.cds) - 3
        removed_indices = mapping[insertion_t0:insertion_t0 + 3]
        if len(removed_indices) != 3 or max(removed_indices) - min(removed_indices) != 2:
            raise DesignError("The stop codon is not contiguous in the genome")
        insertion_boundary0 = transcript_boundary_to_genome0(mapping, insertion_t0, record.strand)
        removed_start0 = min(removed_indices)
        removed_end0 = max(removed_indices) + 1
        removed_gene_sequence = stop
        cds_without_stop = record.cds[:-3]
    elif terminus_upper.startswith("N"):
        if record.cds[:3] != "ATG":
            raise DesignError(f"Transcript CDS does not begin with ATG: {record.cds[:3]}")
        insertion_t0 = cds_start + 3
        insertion_boundary0 = transcript_boundary_to_genome0(mapping, insertion_t0, record.strand)
        removed_start0 = insertion_boundary0
        removed_end0 = insertion_boundary0
        removed_gene_sequence = ""
        cds_without_stop = (
            record.cds[:-3] if record.cds[-3:] in {"TAA", "TAG", "TGA"} else record.cds
        )
    else:
        raise DesignError("Terminus must be N-terminal or C-terminal")

    if record.strand == 1:
        five_start0, five_end0 = insertion_boundary0 - arm_length, insertion_boundary0
        three_start0 = removed_end0
        three_end0 = removed_end0 + arm_length
    else:
        five_start0, five_end0 = insertion_boundary0, insertion_boundary0 + arm_length
        three_start0 = removed_start0 - arm_length
        three_end0 = removed_start0
    if min(five_start0, three_start0) < 0:
        raise DesignError("A homology arm would extend before the chromosome start")

    five_arm = _make_arm(
        name="5-prime homology arm",
        chromosome=record.chromosome,
        start0=five_start0,
        end0=five_end0,
        chromosome_forward_sequence=client.region_sequence(
            record.species, record.chromosome, five_start0, five_end0
        ),
        gene_strand=record.strand,
    )
    three_arm = _make_arm(
        name="3-prime homology arm",
        chromosome=record.chromosome,
        start0=three_start0,
        end0=three_end0,
        chromosome_forward_sequence=client.region_sequence(
            record.species, record.chromosome, three_start0, three_end0
        ),
        gene_strand=record.strand,
    )

    margin = guide_window + 30
    guide_start0 = max(0, insertion_boundary0 - margin)
    guide_end0 = insertion_boundary0 + margin
    guides = enumerate_spcas9_guides(
        chromosome_forward_sequence=client.region_sequence(
            record.species, record.chromosome, guide_start0, guide_end0
        ),
        region_start0=guide_start0,
        insertion_boundary0=insertion_boundary0,
        search_radius=guide_window,
        removed_start0=removed_start0,
        removed_end0=removed_end0,
    )

    return _finalize_result(
        record=record,
        terminus=terminus,
        arm_length=arm_length,
        guide_window=guide_window,
        insertion_boundary0=insertion_boundary0,
        removed_start0=removed_start0,
        removed_end0=removed_end0,
        removed_gene_sequence=removed_gene_sequence,
        cds_without_stop=cds_without_stop,
        edited_cds_without_stop=cds_without_stop,
        five_prime_arm=five_arm,
        three_prime_arm=three_arm,
        guides=guides,
        provenance=[
            f"Gene, transcript, and genomic sequence retrieved live from Ensembl REST ({record.species.assembly}).",
            f"Transcript selected: {record.display_transcript_id}.",
            "Guide/arm rules and SapI adapters follow Bollen et al. 2022 supplementary S1/S3.",
            "Uploaded Addgene #169227 SnapGene sequence parsed and checked against the Bollen supplementary SapI overhang architecture.",
            "Fixed linker-mNeonGreen-stop payload extracted from the uploaded backbone and verified against Bollen supplementary S2.",
        ],
    )


def design_tubb5_fixture(
    *, arm_length: int = 600, guide_window: int = 50
) -> DesignResult:
    if arm_length != 600:
        raise DesignError(
            "The bundled Tubb5 validation fixture is fixed at 600-bp arms; use live Ensembl mode for another length."
        )
    cdna = load_tubb5_cdna()
    if cdna[TUBB5_INSERTION_CDNA0:TUBB5_INSERTION_CDNA0 + 3] != TUBB5_STOP_CODON:
        raise DesignError("Bundled Tubb5 stop-codon coordinate failed validation")
    cds_without_stop = cdna[TUBB5_CDS_START_CDNA0:TUBB5_INSERTION_CDNA0]
    if len(translate(cds_without_stop)) != TUBB5_PROTEIN_LENGTH_AA:
        raise DesignError("Bundled Tubb5 CDS did not translate to 444 aa")

    five_cdna_start0 = TUBB5_INSERTION_CDNA0 - arm_length
    five_gene = cdna[five_cdna_start0:TUBB5_INSERTION_CDNA0]
    three_gene = cdna[
        TUBB5_INSERTION_CDNA0 + 3:TUBB5_INSERTION_CDNA0 + 3 + arm_length
    ]
    five_arm = _make_arm(
        name="5-prime homology arm",
        chromosome=TUBB5_CHROMOSOME,
        start0=TUBB5_INSERTION_BOUNDARY0,
        end0=TUBB5_INSERTION_BOUNDARY0 + arm_length,
        chromosome_forward_sequence=reverse_complement(five_gene),
        gene_strand=TUBB5_STRAND,
    )
    three_arm = _make_arm(
        name="3-prime homology arm",
        chromosome=TUBB5_CHROMOSOME,
        start0=TUBB5_STOP_START0 - arm_length,
        end0=TUBB5_STOP_START0,
        chromosome_forward_sequence=reverse_complement(three_gene),
        gene_strand=TUBB5_STRAND,
    )

    margin = guide_window + 30
    transcript_local_start = TUBB5_INSERTION_CDNA0 - margin
    transcript_local_end = TUBB5_INSERTION_CDNA0 + 3 + margin
    guide_gene_sequence = cdna[transcript_local_start:transcript_local_end]
    guide_region_start0 = TUBB5_INSERTION_BOUNDARY0 + (
        TUBB5_INSERTION_CDNA0 - transcript_local_end
    )
    guides = enumerate_spcas9_guides(
        chromosome_forward_sequence=reverse_complement(guide_gene_sequence),
        region_start0=guide_region_start0,
        insertion_boundary0=TUBB5_INSERTION_BOUNDARY0,
        search_radius=guide_window,
        removed_start0=TUBB5_STOP_START0,
        removed_end0=TUBB5_STOP_END0,
    )
    if not guides:
        raise DesignError("No Tubb5 guide was found in the fixture")
    if guides[0].target_with_pam != "GAGGCAGAAGAGGAGGCCTAAGG":
        raise DesignError("Tubb5 fixture selected an unexpected top guide")

    five_arm, cdna_changes = _tubb5_corrected_five_prime_arm(
        five_arm, cdna=cdna, arm_cdna_start0=five_cdna_start0
    )
    edited_cds_without_stop = _apply_cdna_changes_to_cds(
        cds_without_stop,
        cds_start_cdna0=TUBB5_CDS_START_CDNA0,
        changes=cdna_changes,
    )
    if translate(edited_cds_without_stop) != translate(cds_without_stop):
        raise DesignError("The proposed Tubb5 arm corrections alter the protein sequence")

    # The selected guide's PAM overlaps the endogenous stop codon and is removed
    # by the C-terminal insertion. Following the Bollen guide-selection logic,
    # that destroys the edited-allele target without an additional blocking edit.
    if not guides[0].pam_destroyed or guides[0].blocking_mutation_required:
        raise DesignError("Expected the selected Tubb5 guide PAM to be destroyed by the edit")
    guides[0].blocking_mutation_note = (
        "No extra guide-blocking mutation is required: replacement of the endogenous "
        "TAA stop codon removes one base of the selected guide's AGG PAM."
    )


    record = TranscriptRecord(
        species=TUBB5_SPECIES,
        gene_symbol="Tubb5",
        gene_id=TUBB5_GENE_ID,
        transcript_id=TUBB5_TRANSCRIPT_ID,
        transcript_version=TUBB5_FIXTURE_VERSION,
        chromosome=TUBB5_CHROMOSOME,
        strand=TUBB5_STRAND,
        cdna=cdna,
        cds=cdna[TUBB5_CDS_START_CDNA0:TUBB5_INSERTION_CDNA0 + 3],
        exons=[],
        source="Bundled UCSC/GENCODE mm39 fixture",
        source_release="GENCODE VM snapshot",
    )

    return _finalize_result(
        record=record,
        terminus="C-terminal",
        arm_length=arm_length,
        guide_window=guide_window,
        insertion_boundary0=TUBB5_INSERTION_BOUNDARY0,
        removed_start0=TUBB5_STOP_START0,
        removed_end0=TUBB5_STOP_END0,
        removed_gene_sequence=TUBB5_STOP_CODON,
        cds_without_stop=cds_without_stop,
        edited_cds_without_stop=edited_cds_without_stop,
        five_prime_arm=five_arm,
        three_prime_arm=three_arm,
        guides=guides,
        provenance=[
            "Offline validation sequence: UCSC/GENCODE mm39 ENSMUST00000001566.10 fixture.",
            f"Current Ensembl canonical stable transcript at build time: ENSMUST00000001566.{TUBB5_CURRENT_CANONICAL_VERSION_AT_BUILD}; it remains a 444-aa Tubb5 protein.",
            "The designed coding boundary and 600-bp terminal-exon arms are on the GRCm39/mm39 reference locus.",
            "Guide selection, 14-nt retargeting cutoff, 600-bp arms, and SapI adapters follow Bollen et al. 2022 supplementary S1/S3.",
            "The uploaded Addgene #169227 SnapGene backbone was parsed for full circular in-silico Golden Gate assembly.",
            "The fixed GGGGSAS-mNeonGreen-stop payload was extracted from the backbone and verified against Bollen supplementary S2.",
        ],
        extra_warnings=[
            "Tubb5 has a functionally important C-terminal tail; a C-terminal fluorescent fusion may perturb tubulin interactions or post-translational modification. This test is computational only.",
            "The bundled transcript is ENSMUST00000001566.10, matching the current Ensembl canonical Tubb5-201 listing checked for this build. Live mode should still be rerun before experimental use.",
        ],
    )
