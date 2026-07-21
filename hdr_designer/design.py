from __future__ import annotations

from dataclasses import dataclass, replace
from itertools import combinations

from .backbones import (
    BackboneDefinition,
    SAPI_RECOGNITION_MOTIFS,
    backbone_for_terminus,
    payload_metadata_for,
    synthesis_fragments_for_backbone,
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
from .guides import (
    GUIDE_SAFETY_CUTOFF_NT,
    enumerate_spcas9_guides,
    longest_retained_segment_after_point_mutations,
)
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


SPLICE_EDGE_EXCLUSION_NT = 3
NONCODING_SPLICE_FLANK_EXCLUSION_NT = 6


@dataclass(frozen=True)
class _SynonymousPlan:
    arm_name: str
    codon_start_cdna0: int
    original_codon: str
    altered_codon: str
    amino_acid: str
    # arm index0, genomic position0, cDNA position0, reference base, alternate base
    changes: tuple[tuple[int, int, int, str, str], ...]

    @property
    def change_count(self) -> int:
        return len(self.changes)


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


def _genomic_position0_to_arm_index(
    arm: HomologyArm, genomic_position0: int, strand: int
) -> int | None:
    if not arm.genomic_start0 <= genomic_position0 < arm.genomic_end0:
        return None
    return (
        genomic_position0 - arm.genomic_start0
        if strand == 1
        else arm.genomic_end0 - 1 - genomic_position0
    )


def _synonymous_plans_for_region(
    *,
    arm: HomologyArm,
    strand: int,
    transcript_mapping: list[int],
    cds_start_cdna0: int,
    cds_length: int,
    current_cdna: list[str],
    eligible_genomic_positions0: set[int],
    protected_genomic_positions0: set[int] | None = None,
) -> list[_SynonymousPlan]:
    """Enumerate synonymous codon replacements touching a requested genomic region."""
    protected_genomic_positions0 = protected_genomic_positions0 or set()
    genome_to_cdna = {genomic0: cdna0 for cdna0, genomic0 in enumerate(transcript_mapping)}
    cds_end_cdna0 = cds_start_cdna0 + cds_length
    codon_starts: set[int] = set()
    for genomic0 in eligible_genomic_positions0:
        cdna0 = genome_to_cdna.get(genomic0)
        if cdna0 is None or not cds_start_cdna0 <= cdna0 < cds_end_cdna0:
            continue
        cds_index0 = cdna0 - cds_start_cdna0
        codon_starts.add(cdna0 - (cds_index0 % 3))

    plans: list[_SynonymousPlan] = []
    seen: set[tuple[tuple[int, str], ...]] = set()
    for codon_start_cdna0 in sorted(codon_starts):
        original_codon = "".join(current_cdna[codon_start_cdna0:codon_start_cdna0 + 3])
        if len(original_codon) != 3:
            continue
        amino_acid = CODON_TABLE.get(original_codon)
        # Do not rewrite a terminal stop automatically; stop-context changes require manual design.
        if not amino_acid or amino_acid == "*":
            continue
        for altered_codon, altered_amino_acid in sorted(CODON_TABLE.items()):
            if altered_codon == original_codon or altered_amino_acid != amino_acid:
                continue
            changes: list[tuple[int, int, int, str, str]] = []
            valid = True
            touches_eligible_position = False
            for offset, (reference_base, alternate_base) in enumerate(
                zip(original_codon, altered_codon)
            ):
                if reference_base == alternate_base:
                    continue
                cdna0 = codon_start_cdna0 + offset
                genomic0 = transcript_mapping[cdna0]
                if genomic0 in protected_genomic_positions0:
                    valid = False
                    break
                arm_index0 = _genomic_position0_to_arm_index(arm, genomic0, strand)
                if arm_index0 is None:
                    valid = False
                    break
                if arm.final_gene_oriented_sequence[arm_index0] != reference_base:
                    valid = False
                    break
                touches_eligible_position |= genomic0 in eligible_genomic_positions0
                changes.append(
                    (arm_index0, genomic0, cdna0, reference_base, alternate_base)
                )
            if not valid or not changes or not touches_eligible_position:
                continue
            signature = tuple(sorted((genomic0, alternate) for _, genomic0, _, _, alternate in changes))
            if signature in seen:
                continue
            seen.add(signature)
            plans.append(
                _SynonymousPlan(
                    arm_name=arm.name,
                    codon_start_cdna0=codon_start_cdna0,
                    original_codon=original_codon,
                    altered_codon=altered_codon,
                    amino_acid=amino_acid,
                    changes=tuple(changes),
                )
            )
    return plans


def _splice_edge_exclusion_positions0(record: TranscriptRecord) -> set[int]:
    protected: set[int] = set()
    for exon in record.exons:
        start0 = exon.start1 - 1
        end0 = exon.end1
        protected.update(range(start0, min(end0, start0 + SPLICE_EDGE_EXCLUSION_NT)))
        protected.update(range(max(start0, end0 - SPLICE_EDGE_EXCLUSION_NT), end0))
    return protected


def _noncoding_sapi_exclusion_positions0(record: TranscriptRecord) -> set[int]:
    protected = _splice_edge_exclusion_positions0(record)
    for exon in record.exons:
        start0 = exon.start1 - 1
        end0 = exon.end1
        protected.update(
            range(
                max(0, start0 - NONCODING_SPLICE_FLANK_EXCLUSION_NT),
                start0,
            )
        )
        protected.update(
            range(end0, end0 + NONCODING_SPLICE_FLANK_EXCLUSION_NT)
        )
    return protected


def _longest_homopolymer(sequence: str) -> int:
    longest = 0
    current = 0
    previous = ""
    for base in sequence:
        current = current + 1 if base == previous else 1
        longest = max(longest, current)
        previous = base
    return longest


def _creates_longer_homopolymer(before: str, after: str) -> bool:
    before_longest = _longest_homopolymer(before)
    after_longest = _longest_homopolymer(after)
    return after_longest >= 6 and after_longest > before_longest


def _is_transition(reference: str, alternate: str) -> bool:
    return (reference, alternate) in {
        ("A", "G"),
        ("G", "A"),
        ("C", "T"),
        ("T", "C"),
    }


def _sequence_after_plans(
    arm: HomologyArm, plans: tuple[_SynonymousPlan, ...]
) -> str:
    sequence = list(arm.final_gene_oriented_sequence)
    for plan in plans:
        if plan.arm_name != arm.name:
            continue
        for arm_index0, _, _, _, alternate_base in plan.changes:
            sequence[arm_index0] = alternate_base
    return "".join(sequence)


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
    protein_consequence: str = "",
    pam_before: str = "",
    pam_after: str = "",
    longest_retained_before: int | None = None,
    longest_retained_after: int | None = None,
    automatic: bool = True,
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
        protein_consequence=protein_consequence,
        pam_before=pam_before,
        pam_after=pam_after,
        longest_retained_before=longest_retained_before,
        longest_retained_after=longest_retained_after,
        automatic=automatic,
        reason=reason,
    )
    consequence_note = (
        f"{original_codon}->{altered_codon}, {amino_acid} unchanged."
        if original_codon and altered_codon
        else f"{protein_consequence or 'non-coding sequence'}."
    )
    note = (
        f"{kind}: gene-oriented {reference_base}>{alternate_base} at arm base {index0 + 1} "
        f"(chr{arm.chromosome}:{mutation.genomic_position1:,}); {consequence_note}"
    )
    return replace(
        arm,
        corrected_gene_oriented_sequence=corrected,
        final_sapi_sites=_arm_site_records(corrected),
        mutations=[*arm.mutations, mutation],
        correction_note=" ".join(part for part in [arm.correction_note, note] if part),
    )


def _apply_synonymous_plan(
    arm: HomologyArm,
    plan: _SynonymousPlan,
    *,
    strand: int,
    current_cdna: list[str],
    kind: str,
    reason: str,
    pam_before: str = "",
    pam_after: str = "",
    longest_retained_before: int | None = None,
    longest_retained_after: int | None = None,
) -> HomologyArm:
    if arm.name != plan.arm_name:
        return arm
    for arm_index0, _, cdna0, _, alternate_base in plan.changes:
        arm = _apply_point_mutation(
            arm,
            index0=arm_index0,
            alternate_base=alternate_base,
            strand=strand,
            kind=kind,
            transcript_position1=cdna0 + 1,
            original_codon=plan.original_codon,
            altered_codon=plan.altered_codon,
            amino_acid=f"{plan.amino_acid} unchanged",
            protein_consequence=f"synonymous ({plan.amino_acid})",
            pam_before=pam_before,
            pam_after=pam_after,
            longest_retained_before=longest_retained_before,
            longest_retained_after=longest_retained_after,
            reason=reason,
        )
        current_cdna[cdna0] = alternate_base
    return arm


def _automatic_noncoding_sapi_mutation(
    *,
    arm: HomologyArm,
    motif: str,
    motif_start0: int,
    record: TranscriptRecord,
    transcript_mapping: list[int],
    cds_start_cdna0: int,
    current_cdna: list[str],
) -> HomologyArm | None:
    """Remove a SapI site with one safe, deterministic non-coding substitution."""
    genome_to_cdna = {
        genomic0: cdna0 for cdna0, genomic0 in enumerate(transcript_mapping)
    }
    cds_end_cdna0 = cds_start_cdna0 + len(record.cds)
    protected = _noncoding_sapi_exclusion_positions0(record)
    original_sequence = arm.final_gene_oriented_sequence
    old_site_keys = {
        (str(site["motif"]), int(site["position0"]))
        for site in _arm_site_records(original_sequence)
    }
    candidates: list[
        tuple[tuple[int, int, float, int, int, str], int, int, int | None, str]
    ] = []
    motif_center0 = motif_start0 + (len(motif) - 1) / 2
    for arm_index0 in range(motif_start0, motif_start0 + len(motif)):
        genomic0 = _arm_index_to_genomic_position1(arm, arm_index0, record.strand) - 1
        cdna0 = genome_to_cdna.get(genomic0)
        if cdna0 is not None and cds_start_cdna0 <= cdna0 < cds_end_cdna0:
            continue
        if genomic0 in protected:
            continue
        reference = original_sequence[arm_index0]
        for alternate in "ACGT":
            if alternate == reference:
                continue
            candidate = (
                original_sequence[:arm_index0]
                + alternate
                + original_sequence[arm_index0 + 1:]
            )
            if _creates_longer_homopolymer(original_sequence, candidate):
                continue
            new_site_keys = {
                (str(site["motif"]), int(site["position0"]))
                for site in _arm_site_records(candidate)
            }
            if (motif, motif_start0) in new_site_keys:
                continue
            if not new_site_keys.issubset(old_site_keys):
                continue
            score = (
                0 if cdna0 is None else 1,
                _longest_homopolymer(
                    candidate[max(0, arm_index0 - 6):arm_index0 + 7]
                ),
                abs(arm_index0 - motif_center0),
                0 if _is_transition(reference, alternate) else 1,
                arm_index0,
                alternate,
            )
            candidates.append(
                (score, arm_index0, genomic0, cdna0, alternate)
            )
    if not candidates:
        return None

    _, arm_index0, genomic0, cdna0, alternate = min(
        candidates, key=lambda item: item[0]
    )
    region = "outside the mature transcript" if cdna0 is None else "in an untranslated region"
    arm = _apply_point_mutation(
        arm,
        index0=arm_index0,
        alternate_base=alternate,
        strand=record.strand,
        kind="SapI domestication",
        transcript_position1=cdna0 + 1 if cdna0 is not None else None,
        original_codon="",
        altered_codon="",
        amino_acid="",
        protein_consequence="non-coding (not translated)",
        reason=(
            f"Automatically remove the internal {motif} SapI recognition site with one "
            f"non-coding substitution {region}. Coding bases and the first/last "
            f"{SPLICE_EDGE_EXCLUSION_NT} exonic bases plus "
            f"{NONCODING_SPLICE_FLANK_EXCLUSION_NT} intronic splice-flank bases are "
            "excluded; no new SapI site or "
            "long homopolymer is created. Eligible candidates are ranked by outside-transcript "
            "location, local homopolymer length, distance from the motif center, then transition."
        ),
    )
    if cdna0 is not None:
        current_cdna[cdna0] = alternate
    if arm.mutations[-1].genomic_position1 != genomic0 + 1:
        raise DesignError("Non-coding SapI mutation coordinate validation failed")
    return arm


def _domesticate_sapi_sites(
    *,
    arm: HomologyArm,
    record: TranscriptRecord,
    transcript_mapping: list[int],
    cds_start_cdna0: int,
    current_cdna: list[str],
) -> tuple[HomologyArm, list[str]]:
    """Remove SapI sites with synonymous coding or guarded non-coding changes."""
    warnings: list[str] = []
    unresolved: set[tuple[str, int]] = set()
    while True:
        current_sites = _arm_site_records(arm.final_gene_oriented_sequence)
        pending = [
            site
            for site in current_sites
            if (str(site["motif"]), int(site["position0"])) not in unresolved
        ]
        if not pending:
            break
        site = pending[0]
        motif = str(site["motif"])
        motif_start0 = int(site["position0"])
        eligible_genomic_positions0 = {
            _arm_index_to_genomic_position1(arm, arm_index0, record.strand) - 1
            for arm_index0 in range(motif_start0, motif_start0 + len(motif))
        }
        plans = _synonymous_plans_for_region(
            arm=arm,
            strand=record.strand,
            transcript_mapping=transcript_mapping,
            cds_start_cdna0=cds_start_cdna0,
            cds_length=len(record.cds),
            current_cdna=current_cdna,
            eligible_genomic_positions0=eligible_genomic_positions0,
            protected_genomic_positions0=_splice_edge_exclusion_positions0(record),
        )
        old_site_keys = {
            (str(item["motif"]), int(item["position0"])) for item in current_sites
        }
        original_cds = "".join(
            current_cdna[cds_start_cdna0:cds_start_cdna0 + len(record.cds)]
        )
        viable: list[tuple[tuple[int, int, str], _SynonymousPlan]] = []
        for plan in plans:
            candidate_sequence = _sequence_after_plans(arm, (plan,))
            if _creates_longer_homopolymer(
                arm.final_gene_oriented_sequence, candidate_sequence
            ):
                continue
            new_sites = _arm_site_records(candidate_sequence)
            new_site_keys = {
                (str(item["motif"]), int(item["position0"])) for item in new_sites
            }
            if (motif, motif_start0) in new_site_keys:
                continue
            if not new_site_keys.issubset(old_site_keys):
                continue
            candidate_cdna = list(current_cdna)
            for _, _, cdna0, _, alternate_base in plan.changes:
                candidate_cdna[cdna0] = alternate_base
            candidate_cds = "".join(
                candidate_cdna[cds_start_cdna0:cds_start_cdna0 + len(record.cds)]
            )
            if translate(candidate_cds) != translate(original_cds):
                continue
            changes_outside_site = sum(
                not (motif_start0 <= arm_index0 < motif_start0 + len(motif))
                for arm_index0, _, _, _, _ in plan.changes
            )
            viable.append(
                ((plan.change_count, changes_outside_site, plan.altered_codon), plan)
            )
        if not viable:
            noncoding_arm = _automatic_noncoding_sapi_mutation(
                arm=arm,
                motif=motif,
                motif_start0=motif_start0,
                record=record,
                transcript_mapping=transcript_mapping,
                cds_start_cdna0=cds_start_cdna0,
                current_cdna=current_cdna,
            )
            if noncoding_arm is not None:
                arm = noncoding_arm
                continue
            unresolved.add((motif, motif_start0))
            genomic_positions = sorted(position + 1 for position in eligible_genomic_positions0)
            warnings.append(
                f"Internal SapI site {motif} in {arm.name} at arm base {motif_start0 + 1} "
                f"(genomic bases {genomic_positions[0]}-{genomic_positions[-1]}) could not be "
                "removed by a verified synonymous coding change or an eligible non-coding "
                "single-base substitution; manual intervention is required."
            )
            continue
        plan = min(viable, key=lambda item: item[0])[1]
        arm = _apply_synonymous_plan(
            arm,
            plan,
            strand=record.strand,
            current_cdna=current_cdna,
            kind="SapI domestication",
            reason=(
                f"Automatically remove the internal {motif} SapI recognition site with "
                "the smallest verified synonymous codon replacement; the complete CDS "
                "translation is unchanged and no new SapI site is created."
            ),
        )
    return arm, warnings


def _donor_target_after_edits(
    guide: GuideCandidate,
    *,
    gene_strand: int,
    arms: tuple[HomologyArm, HomologyArm],
    additional_plans: tuple[_SynonymousPlan, ...] = (),
) -> str:
    """Apply gene-oriented arm substitutions to a guide-oriented target sequence."""
    target = list(guide.target_with_pam)
    edits: dict[int, str] = {}
    for arm in arms:
        for mutation in arm.mutations:
            edits[mutation.genomic_position1 - 1] = mutation.alternate_base
    for plan in additional_plans:
        for _, genomic0, _, _, alternate_base in plan.changes:
            edits[genomic0] = alternate_base

    for genomic0, alternate_gene_base in edits.items():
        if not guide.target_start0 <= genomic0 < guide.target_end0:
            continue
        target_index0 = (
            genomic0 - guide.target_start0
            if guide.chromosome_strand == "+"
            else guide.target_end0 - 1 - genomic0
        )
        alternate_chromosome_base = (
            alternate_gene_base
            if gene_strand == 1
            else reverse_complement(alternate_gene_base)
        )
        target[target_index0] = (
            alternate_chromosome_base
            if guide.chromosome_strand == "+"
            else reverse_complement(alternate_chromosome_base)
        )
    return "".join(target)


def _mutated_genomic_positions0(
    arms: tuple[HomologyArm, HomologyArm],
    additional_plans: tuple[_SynonymousPlan, ...] = (),
) -> list[int]:
    positions = {
        mutation.genomic_position1 - 1
        for arm in arms
        for mutation in arm.mutations
    }
    positions.update(
        genomic0
        for plan in additional_plans
        for _, genomic0, _, _, _ in plan.changes
    )
    return sorted(positions)


def _retained_after_donor_edits(
    guide: GuideCandidate,
    *,
    insertion_boundary0: int,
    removed_start0: int,
    removed_end0: int,
    arms: tuple[HomologyArm, HomologyArm],
    additional_plans: tuple[_SynonymousPlan, ...] = (),
) -> int:
    return longest_retained_segment_after_point_mutations(
        target_start0=guide.target_start0,
        target_end0=guide.target_end0,
        insertion_boundary0=insertion_boundary0,
        removed_start0=removed_start0,
        removed_end0=removed_end0,
        mutated_genomic_positions0=_mutated_genomic_positions0(arms, additional_plans),
    )


def _guide_plan_preserves_sapi_sites(
    arms: tuple[HomologyArm, HomologyArm], plans: tuple[_SynonymousPlan, ...]
) -> bool:
    for arm in arms:
        before = {
            (str(site["motif"]), int(site["position0"]))
            for site in _arm_site_records(arm.final_gene_oriented_sequence)
        }
        after = {
            (str(site["motif"]), int(site["position0"]))
            for site in _arm_site_records(_sequence_after_plans(arm, plans))
        }
        if not after.issubset(before):
            return False
    return True


def _design_guide_blocking_mutations(
    *,
    guide: GuideCandidate,
    five_prime_arm: HomologyArm,
    three_prime_arm: HomologyArm,
    record: TranscriptRecord,
    transcript_mapping: list[int],
    cds_start_cdna0: int,
    current_cdna: list[str],
    insertion_boundary0: int,
    removed_start0: int,
    removed_end0: int,
) -> tuple[HomologyArm, HomologyArm, list[str]]:
    """Protect the edited allele with verified synonymous PAM/seed changes."""
    warnings: list[str] = []
    arms = (five_prime_arm, three_prime_arm)
    current_target = _donor_target_after_edits(
        guide, gene_strand=record.strand, arms=arms
    )
    current_retained = _retained_after_donor_edits(
        guide,
        insertion_boundary0=insertion_boundary0,
        removed_start0=removed_start0,
        removed_end0=removed_end0,
        arms=arms,
    )
    current_pam = current_target[-3:]
    current_pam_functional = current_target[21:23] == "GG"
    guide.final_pam = current_pam
    guide.final_pam_destroyed = guide.pam_destroyed or not current_pam_functional
    guide.final_longest_retained_segment = current_retained
    guide.blocking_mutation_required = (
        not guide.final_pam_destroyed
        and current_retained > GUIDE_SAFETY_CUTOFF_NT
    )

    if not guide.blocking_mutation_required:
        if guide.pam_destroyed:
            return five_prime_arm, three_prime_arm, warnings
        if current_target != guide.target_with_pam:
            guide.blocking_mutation_note = (
                "No additional guide-blocking mutation is required: an automatically "
                f"designed donor correction changes the PAM from {guide.pam} to {current_pam} "
                f"or reduces the longest retained segment from {guide.longest_retained_segment} "
                f"to {current_retained} nt."
            )
        return five_prime_arm, three_prime_arm, warnings

    eligible_genomic_positions0 = set(range(guide.target_start0, guide.target_end0))
    existing_mutated_positions0 = set(_mutated_genomic_positions0(arms))
    plans: list[_SynonymousPlan] = []
    for arm in arms:
        plans.extend(
            _synonymous_plans_for_region(
                arm=arm,
                strand=record.strand,
                transcript_mapping=transcript_mapping,
                cds_start_cdna0=cds_start_cdna0,
                cds_length=len(record.cds),
                current_cdna=current_cdna,
                eligible_genomic_positions0=eligible_genomic_positions0,
                protected_genomic_positions0=_splice_edge_exclusion_positions0(record),
            )
        )
    plans = [
        plan
        for plan in plans
        if not existing_mutated_positions0.intersection(
            genomic0 for _, genomic0, _, _, _ in plan.changes
        )
    ]

    original_cds = "".join(
        current_cdna[cds_start_cdna0:cds_start_cdna0 + len(record.cds)]
    )
    viable: list[
        tuple[tuple[int, int, int, int, tuple[tuple[int, str], ...]], tuple[_SynonymousPlan, ...], str, int]
    ] = []
    for plan_count in range(1, min(3, len(plans)) + 1):
        for selected in combinations(plans, plan_count):
            codon_keys = {(plan.arm_name, plan.codon_start_cdna0) for plan in selected}
            if len(codon_keys) != len(selected):
                continue
            merged: dict[int, str] = {}
            conflict = False
            for plan in selected:
                for _, genomic0, _, _, alternate_base in plan.changes:
                    if genomic0 in merged and merged[genomic0] != alternate_base:
                        conflict = True
                        break
                    merged[genomic0] = alternate_base
                if conflict:
                    break
            if conflict or not _guide_plan_preserves_sapi_sites(arms, selected):
                continue
            if any(
                _creates_longer_homopolymer(
                    arm.final_gene_oriented_sequence,
                    _sequence_after_plans(arm, selected),
                )
                for arm in arms
            ):
                continue
            candidate_cdna = list(current_cdna)
            for plan in selected:
                for _, _, cdna0, _, alternate_base in plan.changes:
                    candidate_cdna[cdna0] = alternate_base
            candidate_cds = "".join(
                candidate_cdna[cds_start_cdna0:cds_start_cdna0 + len(record.cds)]
            )
            if translate(candidate_cds) != translate(original_cds):
                continue
            target_after = _donor_target_after_edits(
                guide,
                gene_strand=record.strand,
                arms=arms,
                additional_plans=selected,
            )
            pam_after = target_after[-3:]
            pam_destroyed = target_after[21:23] != "GG"
            retained_after = _retained_after_donor_edits(
                guide,
                insertion_boundary0=insertion_boundary0,
                removed_start0=removed_start0,
                removed_end0=removed_end0,
                arms=arms,
                additional_plans=selected,
            )
            if not pam_destroyed and retained_after > GUIDE_SAFETY_CUTOFF_NT:
                continue
            target_indexes = [
                (
                    genomic0 - guide.target_start0
                    if guide.chromosome_strand == "+"
                    else guide.target_end0 - 1 - genomic0
                )
                for genomic0 in merged
            ]
            closest_to_pam = min(
                (0 if index0 >= 20 else 20 - index0) for index0 in target_indexes
            )
            signature = tuple(sorted(merged.items()))
            score = (
                0 if pam_destroyed else 1,
                len(merged),
                closest_to_pam,
                retained_after,
                signature,
            )
            viable.append((score, selected, pam_after, retained_after))

    if not viable:
        guide.blocking_mutation_note = (
            "A guide-blocking mutation is still required. No synonymous coding change "
            "within the retained target could destroy the PAM or satisfy the 14-nt cutoff; "
            "noncoding changes are not released automatically. Manual sequence design or "
            "another guide is required, so no sequence-complete output is released."
        )
        warnings.append(guide.blocking_mutation_note)
        return five_prime_arm, three_prime_arm, warnings

    _, selected_plans, final_pam, final_retained = min(viable, key=lambda item: item[0])
    pam_destroyed = final_pam[1:3] != "GG"
    reason = (
        f"Automatically protect the edited allele by changing the selected guide PAM from "
        f"{current_pam} to {final_pam} with a verified synonymous codon replacement."
        if pam_destroyed
        else (
            "Automatically protect the edited allele with the smallest verified synonymous "
            f"seed change set, reducing the longest retained target segment from "
            f"{current_retained} to {final_retained} nt (cutoff <= {GUIDE_SAFETY_CUTOFF_NT} nt)."
        )
    )
    for plan in selected_plans:
        if plan.arm_name == five_prime_arm.name:
            five_prime_arm = _apply_synonymous_plan(
                five_prime_arm,
                plan,
                strand=record.strand,
                current_cdna=current_cdna,
                kind="Guide blocking",
                reason=reason,
                pam_before=current_pam,
                pam_after=final_pam,
                longest_retained_before=current_retained,
                longest_retained_after=final_retained,
            )
        elif plan.arm_name == three_prime_arm.name:
            three_prime_arm = _apply_synonymous_plan(
                three_prime_arm,
                plan,
                strand=record.strand,
                current_cdna=current_cdna,
                kind="Guide blocking",
                reason=reason,
                pam_before=current_pam,
                pam_after=final_pam,
                longest_retained_before=current_retained,
                longest_retained_after=final_retained,
            )

    final_arms = (five_prime_arm, three_prime_arm)
    verified_target = _donor_target_after_edits(
        guide, gene_strand=record.strand, arms=final_arms
    )
    verified_retained = _retained_after_donor_edits(
        guide,
        insertion_boundary0=insertion_boundary0,
        removed_start0=removed_start0,
        removed_end0=removed_end0,
        arms=final_arms,
    )
    if verified_target[-3:] != final_pam or verified_retained != final_retained:
        raise DesignError("Guide-blocking mutation verification did not reproduce its prediction")
    guide.final_pam = final_pam
    guide.final_pam_destroyed = guide.pam_destroyed or verified_target[21:23] != "GG"
    guide.final_longest_retained_segment = verified_retained
    guide.blocking_mutation_required = (
        not guide.final_pam_destroyed
        and verified_retained > GUIDE_SAFETY_CUTOFF_NT
    )
    if guide.blocking_mutation_required:
        raise DesignError("Automatically designed guide-blocking mutation did not pass the safety gate")
    guide.blocking_mutation_note = (
        f"Automatic synonymous guide blocking applied: PAM {current_pam}->{final_pam}; "
        f"longest retained segment {current_retained}->{verified_retained} nt. "
        "The complete CDS translation is unchanged and no SapI site was created."
    )
    return five_prime_arm, three_prime_arm, warnings


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
    backbone_definition: BackboneDefinition | None = None,
    extra_warnings: list[str] | None = None,
) -> DesignResult:
    is_c_terminal = terminus.upper().startswith("C")
    definition = backbone_definition or backbone_for_terminus(terminus)
    if definition.terminus != ("C-terminal" if is_c_terminal else "N-terminal"):
        raise DesignError(
            f"The selected {definition.terminus} backbone cannot be used for {terminus} tagging"
        )
    top = guides[0] if guides else None
    final_sapi_count = len(five_prime_arm.final_sapi_sites) + len(three_prime_arm.final_sapi_sites)

    donor_payload: dict[str, object] = payload_metadata_for(definition)
    cloning_fragments: dict[str, object] = {}
    primer_tails: dict[str, str] = {}
    edited_cds = ""
    fusion_protein = ""
    junctions: dict[str, str] = {}

    guide_safe = bool(top and not top.blocking_mutation_required)
    can_release = bool(top and guide_safe and final_sapi_count == 0)
    if can_release and top:
        cloning_fragments = synthesis_fragments_for_backbone(
            definition,
            target_with_pam=top.target_with_pam,
            uha=five_prime_arm.final_gene_oriented_sequence,
            dha=three_prime_arm.final_gene_oriented_sequence,
        )
        primer_tails = {
            "UHA_forward_5prime_tail": (
                definition.uha_forward_primer_tail_prefix + top.target_with_pam
            ),
            "UHA_reverse_5prime_tail": definition.uha_reverse_primer_tail,
            "DHA_forward_5prime_tail": definition.dha_forward_primer_tail,
            "DHA_reverse_5prime_tail": (
                definition.dha_reverse_primer_tail_prefix
                + reverse_complement(top.target_with_pam)
            ),
            "note": (
                "Append each tail to a separately designed locus-specific annealing sequence; "
                "this prototype does not design PCR annealing regions."
            ),
        }
        payload_coding = str(donor_payload["payload_coding_sequence"])
        edited_cds = (
            edited_cds_without_stop + payload_coding
            if is_c_terminal
            else edited_cds_without_stop[:3] + payload_coding + edited_cds_without_stop[3:]
        )
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
    backbone_info = cloning_fragments.get("uploaded_backbone", {})
    backbone_verified = bool(
        backbone_info
        and backbone_info.get("addgene_id") == definition.addgene_id
        and backbone_info.get("length_nt") == definition.expected_length_nt
        and backbone_info.get("topology") == "circular"
        and backbone_info.get("sapi_site_count") == 4
        and backbone_info.get("payload_sequence_verified") is True
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
        cloning_fragments.get("assembled_plasmid_length_nt")
        and cloning_fragments.get("assembled_plasmid_topology") == "circular"
        and cloning_fragments.get("assembled_plasmid_sapi_site_count") == 0
        and junctions_verified
    )

    validations: list[dict[str, str]] = [
        {
            "check": "Reference coding frame",
            "status": "PASS" if len(cds_without_stop) % 3 == 0 else "FAIL",
            "detail": f"Native coding sequence without the terminal stop is {len(cds_without_stop)} nt ({len(cds_without_stop) // 3} aa).",
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
                    "status": "PASS" if top.target_destroyed else "INFO",
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
                "status": "PASS" if can_release else "BLOCKED",
                "detail": (
                    f"Exact supplementary S1 {definition.terminus} adapters and "
                    f"{'/'.join(definition.overhangs.values())} overhangs were applied."
                    if can_release
                    else "Fragments are withheld until guide-retargeting and internal-SapI gates pass."
                ),
            },
            {
                "check": f"Uploaded Addgene #{definition.addgene_id} backbone",
                "status": "PASS" if backbone_verified else "BLOCKED",
                "detail": (
                    f"Parsed {backbone_info.get('snapgene_file')} ({backbone_info.get('length_nt')} bp, "
                    f"{backbone_info.get('topology')}); four SapI sites yield "
                    f"{'/'.join(definition.overhangs.values())}, and the extracted "
                    f"{definition.payload_length_nt}-bp payload matches the verified sequence."
                    if backbone_verified
                    else "The uploaded fixed-backbone sequence did not pass structural verification."
                ),
            },
            {
                "check": "Full circular plasmid assembly",
                "status": "PASS" if plasmid_assembly_verified else "BLOCKED",
                "detail": (
                    f"Reconstructed a {cloning_fragments.get('assembled_plasmid_length_nt')}-bp circular plasmid; "
                    "all four ligation junctions match and no SapI recognition site remains."
                    if plasmid_assembly_verified
                    else "No fully verified circular Golden Gate product was released."
                ),
            },
            {
                "check": "Fusion translation",
                "status": "PASS" if fusion_protein else "BLOCKED",
                "detail": (
                    f"Predicted fusion is {len(fusion_protein)} aa; linker {donor_payload.get('linker_peptide', '')}; "
                    f"{donor_payload.get('tag_name', 'tag')} {donor_payload.get('tag_length_aa', '')} aa; linker-plus-tag payload "
                    f"{donor_payload.get('payload_peptide_length_aa', '')} aa."
                    if fusion_protein
                    else "No final fusion translation released."
                ),
            },
        ]
    )

    sequence_complete = bool(
        can_release and fusion_protein and backbone_verified and plasmid_assembly_verified
    )
    status = (
        "SEQUENCE-COMPLETE COMPUTATIONAL DESIGN"
        if sequence_complete
        else "DESIGN BLOCKED - NO SEQUENCE-COMPLETE OUTPUT"
    )
    warnings = list(extra_warnings or [])
    warnings.extend(
        [
            "No genome-wide off-target analysis is performed, as requested.",
            "Reference-genome sequence only; strain, cell-line, and clone-specific variants are not assessed.",
            "No experimental activity score is calculated. After the two primary Bollen priorities, ranking uses only GC/poly-T heuristics.",
            "Sequence-complete means the internal computational gates and uploaded-backbone assembly simulation passed; independently verify every sequence and plasmid junction before experimental use.",
            "Custom backbone uploads are accepted only when all four SapI sites, the supported overhang order, the GGGGSAS linker, and the payload reading frame validate.",
        ]
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
        backbone_name=definition.name,
        backbone_addgene_id=definition.addgene_id,
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
            "Rank: nick distance first, then target disruption and basic spacer heuristics. "
            "A nearer guide remains preferred when it can be protected with a verified "
            "synonymous donor mutation; no quantitative on-target or off-target model is run."
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
        custom_backbones_supported=True,
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
    backbone_definition: BackboneDefinition | None = None,
) -> DesignResult:
    if species_key not in SPECIES:
        raise DesignError(f"Unsupported species: {species_key}")
    if arm_length < 100:
        raise DesignError("Homology arms shorter than 100 bp are not supported")
    terminus_upper = terminus.upper()
    if not (terminus_upper.startswith("C") or terminus_upper.startswith("N")):
        raise DesignError("Terminus must be N-terminal or C-terminal")
    selected_backbone = backbone_definition or backbone_for_terminus(terminus)
    expected_terminus = "C-terminal" if terminus_upper.startswith("C") else "N-terminal"
    if selected_backbone.terminus != expected_terminus:
        raise DesignError(
            f"The selected {selected_backbone.terminus} backbone cannot be used for {terminus} tagging"
        )
    client = client or EnsemblClient()
    record = client.transcript_record(SPECIES[species_key], gene, transcript_id)
    backbone_source = (
        f"Custom SnapGene backbone {selected_backbone.dna_path.name}"
        if selected_backbone.is_custom
        else f"Uploaded Addgene #{selected_backbone.addgene_id} SnapGene sequence"
    )
    payload_source = (
        "Custom payload extracted between the verified inner SapI cuts"
        if selected_backbone.is_custom
        else f"Fixed {selected_backbone.terminus} mNeonGreen payload extracted from the uploaded backbone"
    )
    mapping = build_transcript_genome_map(record)
    cds_start = _find_unique_cds(record.cdna, record.cds)
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

    current_cdna = list(record.cdna)
    mutation_warnings: list[str] = []
    five_arm, warnings = _domesticate_sapi_sites(
        arm=five_arm,
        record=record,
        transcript_mapping=mapping,
        cds_start_cdna0=cds_start,
        current_cdna=current_cdna,
    )
    mutation_warnings.extend(warnings)
    three_arm, warnings = _domesticate_sapi_sites(
        arm=three_arm,
        record=record,
        transcript_mapping=mapping,
        cds_start_cdna0=cds_start,
        current_cdna=current_cdna,
    )
    mutation_warnings.extend(warnings)

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
    if guides:
        five_arm, three_arm, warnings = _design_guide_blocking_mutations(
            guide=guides[0],
            five_prime_arm=five_arm,
            three_prime_arm=three_arm,
            record=record,
            transcript_mapping=mapping,
            cds_start_cdna0=cds_start,
            current_cdna=current_cdna,
            insertion_boundary0=insertion_boundary0,
            removed_start0=removed_start0,
            removed_end0=removed_end0,
        )
        mutation_warnings.extend(warnings)

    edited_cds_without_stop = "".join(
        current_cdna[cds_start:cds_start + len(cds_without_stop)]
    )
    if translate(edited_cds_without_stop) != translate(cds_without_stop):
        raise DesignError("Automatic donor mutations changed the reference protein sequence")

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
        edited_cds_without_stop=edited_cds_without_stop,
        five_prime_arm=five_arm,
        three_prime_arm=three_arm,
        guides=guides,
        provenance=[
            f"Gene, transcript, and genomic sequence retrieved live from Ensembl REST ({record.species.assembly}).",
            f"Transcript selected: {record.display_transcript_id}.",
            "Guide/arm rules and SapI adapters follow Bollen et al. 2022 supplementary S1/S3.",
            f"{backbone_source} parsed and checked against the Bollen supplementary {selected_backbone.terminus} SapI overhang architecture.",
            f"{payload_source} and verified by length, reading frame, and SHA-256 digest.",
            "Coding SapI sites and selected-guide retargeting are evaluated with generic synonymous-codon search and full-CDS translation checks.",
        ],
        backbone_definition=selected_backbone,
        extra_warnings=mutation_warnings,
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
