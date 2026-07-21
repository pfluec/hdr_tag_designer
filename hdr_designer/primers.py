from __future__ import annotations

from dataclasses import dataclass
from hashlib import sha256
from typing import Any

import primer3

from .sequence import gc_percent, reverse_complement


PRIMER_RULESET = (
    f"Primer3 2.6.1 via primer3-py {primer3.__version__}; "
    "HDR genotyping constraints v1"
)
PRIMER_MIN_SIZE = 18
PRIMER_OPT_SIZE = 20
PRIMER_MAX_SIZE = 27
PRIMER_MIN_TM = 57.0
PRIMER_OPT_TM = 60.0
PRIMER_MAX_TM = 63.0
PRIMER_MIN_GC = 35.0
PRIMER_MAX_GC = 65.0
PRIMER_MAX_TM_DIFFERENCE = 3.0
PRIMER_MAX_STRUCTURE_TM = 45.0
PAYLOAD_JUNCTION_STANDOFF_NT = 150


def _endpoint_primer(sequence: str, orientation: str) -> tuple[_PrimerCandidate, list[str]]:
    """Choose an annealing sequence anchored exactly at one end of an arm."""
    region_start0 = 0 if orientation == "forward" else max(0, len(sequence) - 60)
    region_end0 = min(len(sequence), 60) if orientation == "forward" else len(sequence)
    candidates = _candidate_primers(
        sequence,
        region_start0=region_start0,
        region_end0=region_end0,
        orientation=orientation,
        preferred_boundary0=0 if orientation == "forward" else len(sequence),
    )
    anchored = [
        candidate
        for candidate in candidates
        if (candidate.bind_start0 == 0 if orientation == "forward" else candidate.bind_end0 == len(sequence))
    ]
    if anchored:
        return anchored[0], []

    length = min(PRIMER_OPT_SIZE, len(sequence))
    start0 = 0 if orientation == "forward" else len(sequence) - length
    binding = sequence[start0:start0 + length]
    primer_sequence = binding if orientation == "forward" else reverse_complement(binding)
    tm = float(primer3.calc_tm(primer_sequence))
    hairpin_tm = float(primer3.calc_hairpin(primer_sequence).tm)
    homodimer_tm = float(primer3.calc_homodimer(primer_sequence).tm)
    fallback = _PrimerCandidate(
        sequence=primer_sequence,
        bind_start0=start0,
        bind_end0=start0 + length,
        orientation=orientation,
        tm_c=round(tm, 2),
        gc_percent=gc_percent(primer_sequence),
        hairpin_tm_c=round(hairpin_tm, 2),
        homodimer_tm_c=round(homodimer_tm, 2),
        score=999.0,
    )
    return fallback, [
        "No arm-end annealing sequence passed every default Primer3 constraint; "
        "the displayed 20-nt endpoint fallback requires manual review."
    ]


def design_homology_arm_cloning_primers(
    *,
    five_prime_arm: Any,
    three_prime_arm: Any,
    tails: dict[str, str],
) -> dict[str, Any]:
    """Complete the four Bollen cloning-primer tails with arm-end annealing regions."""
    arm_specs = (
        ("UHA", five_prime_arm, "UHA_forward_5prime_tail", "UHA_reverse_5prime_tail"),
        ("DHA", three_prime_arm, "DHA_forward_5prime_tail", "DHA_reverse_5prime_tail"),
    )
    result: dict[str, Any] = {
        "status": "PASS",
        "ruleset": PRIMER_RULESET,
        "primers": {},
        "warnings": [],
    }
    for arm_key, arm, forward_tail_key, reverse_tail_key in arm_specs:
        final_sequence = arm.final_gene_oriented_sequence
        forward, forward_warnings = _endpoint_primer(final_sequence, "forward")
        reverse, reverse_warnings = _endpoint_primer(final_sequence, "reverse")
        pair = (("forward", forward, tails[forward_tail_key]), ("reverse", reverse, tails[reverse_tail_key]))
        covered_positions1 = set(range(1, forward.bind_end0 + 1)) | set(
            range(reverse.bind_start0 + 1, len(final_sequence) + 1)
        )
        internal_mutations = [
            mutation
            for mutation in arm.mutations
            if mutation.arm_position1 not in covered_positions1
        ]
        for role, candidate, tail in pair:
            full_sequence = tail + candidate.sequence
            candidate_positions1 = set(
                range(candidate.bind_start0 + 1, candidate.bind_end0 + 1)
            )
            result["primers"][f"{arm_key}_{role}"] = {
                "name": f"{arm_key}_{role}_cloning_primer",
                "arm": arm_key,
                "orientation": role,
                "tail_sequence_5to3": tail,
                "annealing_sequence_5to3": candidate.sequence,
                "full_sequence_5to3": full_sequence,
                "annealing_length_nt": len(candidate.sequence),
                "full_length_nt": len(full_sequence),
                "annealing_tm_c": candidate.tm_c,
                "annealing_gc_percent": candidate.gc_percent,
                "annealing_hairpin_tm_c": candidate.hairpin_tm_c,
                "annealing_homodimer_tm_c": candidate.homodimer_tm_c,
                "arm_binding_interval_1based": f"{candidate.bind_start0 + 1}-{candidate.bind_end0}",
                "incorporated_mutations": [
                    mutation.kind
                    for mutation in arm.mutations
                    if mutation.arm_position1 in candidate_positions1
                ],
            }
        warnings = [f"{arm_key} forward: {item}" for item in forward_warnings]
        warnings.extend(f"{arm_key} reverse: {item}" for item in reverse_warnings)
        if internal_mutations:
            warnings.append(
                f"{arm_key} contains {len(internal_mutations)} required internal mutation(s) "
                "outside the primer annealing regions. Endpoint PCR primers alone will not "
                "introduce them; use a synthesized arm or an additional mutagenesis strategy."
            )
            result[f"{arm_key.lower()}_internal_mutations_not_encoded_by_primers"] = [
                {
                    "kind": mutation.kind,
                    "arm_position1": mutation.arm_position1,
                    "change": f"{mutation.reference_base}>{mutation.alternate_base}",
                }
                for mutation in internal_mutations
            ]
        if warnings:
            result["status"] = "WARNING"
            result["warnings"].extend(warnings)
    return result


@dataclass(frozen=True)
class _PrimerCandidate:
    sequence: str
    bind_start0: int
    bind_end0: int
    orientation: str
    tm_c: float
    gc_percent: float
    hairpin_tm_c: float
    homodimer_tm_c: float
    score: float


def _longest_homopolymer(sequence: str) -> int:
    longest = current = 0
    previous = ""
    for base in sequence:
        current = current + 1 if base == previous else 1
        longest = max(longest, current)
        previous = base
    return longest


def _candidate_primers(
    template: str,
    *,
    region_start0: int,
    region_end0: int,
    orientation: str,
    preferred_boundary0: int,
    limit: int = 120,
) -> list[_PrimerCandidate]:
    candidates: list[_PrimerCandidate] = []
    for length in range(PRIMER_MIN_SIZE, PRIMER_MAX_SIZE + 1):
        for start0 in range(region_start0, region_end0 - length + 1):
            binding = template[start0:start0 + length]
            if set(binding) != set(binding).intersection("ACGT"):
                continue
            sequence = binding if orientation == "forward" else reverse_complement(binding)
            gc = gc_percent(sequence)
            if not PRIMER_MIN_GC <= gc <= PRIMER_MAX_GC:
                continue
            if _longest_homopolymer(sequence) > 4:
                continue
            tm = float(primer3.calc_tm(sequence))
            if not PRIMER_MIN_TM <= tm <= PRIMER_MAX_TM:
                continue
            hairpin_tm = float(primer3.calc_hairpin(sequence).tm)
            homodimer_tm = float(primer3.calc_homodimer(sequence).tm)
            if hairpin_tm >= PRIMER_MAX_STRUCTURE_TM or homodimer_tm >= PRIMER_MAX_STRUCTURE_TM:
                continue
            three_prime_gc = sequence[-5:].count("G") + sequence[-5:].count("C")
            if three_prime_gc > 4:
                continue
            boundary_distance = (
                abs(start0 - preferred_boundary0)
                if orientation == "forward"
                else abs(start0 + length - preferred_boundary0)
            )
            score = (
                abs(tm - PRIMER_OPT_TM) * 3
                + abs(length - PRIMER_OPT_SIZE)
                + boundary_distance * 0.03
                + max(0.0, hairpin_tm - 35.0)
                + max(0.0, homodimer_tm - 35.0)
            )
            candidates.append(
                _PrimerCandidate(
                    sequence=sequence,
                    bind_start0=start0,
                    bind_end0=start0 + length,
                    orientation=orientation,
                    tm_c=round(tm, 2),
                    gc_percent=gc,
                    hairpin_tm_c=round(hairpin_tm, 2),
                    homodimer_tm_c=round(homodimer_tm, 2),
                    score=score,
                )
            )
    return sorted(candidates, key=lambda item: (item.score, item.sequence))[:limit]


def _junction_forward_candidates(
    template: str,
    *,
    junction0: int,
) -> list[_PrimerCandidate]:
    candidates: list[_PrimerCandidate] = []
    for length in range(PRIMER_MIN_SIZE, PRIMER_MAX_SIZE + 1):
        for left_bases in range(6, length - 9):
            start0 = junction0 - left_bases
            end0 = start0 + length
            if start0 < 0 or end0 > len(template) or not start0 < junction0 < end0:
                continue
            binding = template[start0:end0]
            sequence = binding
            gc = gc_percent(sequence)
            if not PRIMER_MIN_GC <= gc <= PRIMER_MAX_GC:
                continue
            if _longest_homopolymer(sequence) > 4:
                continue
            tm = float(primer3.calc_tm(sequence))
            if not PRIMER_MIN_TM <= tm <= PRIMER_MAX_TM:
                continue
            hairpin_tm = float(primer3.calc_hairpin(sequence).tm)
            homodimer_tm = float(primer3.calc_homodimer(sequence).tm)
            if hairpin_tm >= PRIMER_MAX_STRUCTURE_TM or homodimer_tm >= PRIMER_MAX_STRUCTURE_TM:
                continue
            score = (
                abs(tm - PRIMER_OPT_TM) * 3
                + abs(length - PRIMER_OPT_SIZE)
                + abs(left_bases - 8) * 0.5
                + max(0.0, hairpin_tm - 35.0)
                + max(0.0, homodimer_tm - 35.0)
            )
            candidates.append(
                _PrimerCandidate(
                    sequence=sequence,
                    bind_start0=start0,
                    bind_end0=end0,
                    orientation="forward",
                    tm_c=round(tm, 2),
                    gc_percent=gc,
                    hairpin_tm_c=round(hairpin_tm, 2),
                    homodimer_tm_c=round(homodimer_tm, 2),
                    score=score,
                )
            )
    return sorted(candidates, key=lambda item: (item.score, item.sequence))[:120]


def _occurrence_count(sequence: str, query: str) -> int:
    return sequence.count(query) + sequence.count(reverse_complement(query))


def _select_pair(
    *,
    assay: str,
    template: str,
    left_candidates: list[_PrimerCandidate],
    right_candidates: list[_PrimerCandidate],
    assembled_plasmid: str,
    genomic_sides: set[str],
    target_product_size: int,
    priority_side: str | None = None,
) -> tuple[_PrimerCandidate, _PrimerCandidate, dict[str, Any]] | None:
    viable: list[
        tuple[int, float, _PrimerCandidate, _PrimerCandidate, dict[str, Any]]
    ] = []
    for left_rank, left in enumerate(left_candidates):
        if "left" in genomic_sides and _occurrence_count(assembled_plasmid, left.sequence):
            continue
        if _occurrence_count(template, left.sequence) != 1:
            continue
        for right_rank, right in enumerate(right_candidates):
            if right.bind_end0 <= left.bind_start0:
                continue
            if "right" in genomic_sides and _occurrence_count(assembled_plasmid, right.sequence):
                continue
            if _occurrence_count(template, right.sequence) != 1:
                continue
            tm_difference = abs(left.tm_c - right.tm_c)
            if tm_difference > PRIMER_MAX_TM_DIFFERENCE:
                continue
            heterodimer_tm = float(
                primer3.calc_heterodimer(left.sequence, right.sequence).tm
            )
            if heterodimer_tm >= PRIMER_MAX_STRUCTURE_TM:
                continue
            product_size = right.bind_end0 - left.bind_start0
            score = (
                left.score
                + right.score
                + tm_difference * 4
                + abs(product_size - target_product_size) * 0.01
                + max(0.0, heterodimer_tm - 35.0)
            )
            viable.append(
                (
                    (
                        left_rank
                        if priority_side == "left"
                        else right_rank if priority_side == "right" else 0
                    ),
                    score,
                    left,
                    right,
                    {
                        "assay": assay,
                        "status": "PASS",
                        "product_size_bp": product_size,
                        "amplicon_sequence_5to3": template[
                            left.bind_start0:right.bind_end0
                        ],
                        "tm_difference_c": round(tm_difference, 2),
                        "heterodimer_tm_c": round(heterodimer_tm, 2),
                    },
                )
            )
    if not viable:
        return None
    _, _, left, right, metadata = min(viable, key=lambda item: (item[0], item[1]))
    return left, right, metadata


def _genomic_interval(
    *,
    region_start0: int,
    region_end0: int,
    offset_start0: int,
    offset_end0: int,
    gene_strand: int,
) -> str:
    if gene_strand == 1:
        start0 = region_start0 + offset_start0
        end0 = region_start0 + offset_end0
    else:
        start0 = region_end0 - offset_end0
        end0 = region_end0 - offset_start0
    return f"{start0 + 1}-{end0}"


def _primer_dict(
    candidate: _PrimerCandidate,
    *,
    name: str,
    source: str,
    source_offset0: int,
    region_interval0: tuple[int, int] | None,
    gene_strand: int,
    payload_length: int,
) -> dict[str, Any]:
    relative_start0 = candidate.bind_start0 - source_offset0
    relative_end0 = candidate.bind_end0 - source_offset0
    result: dict[str, Any] = {
        "name": name,
        "sequence_5to3": candidate.sequence,
        "orientation": candidate.orientation,
        "reference_sequence_strand": (
            "+"
            if (gene_strand == 1 and candidate.orientation == "forward")
            or (gene_strand == -1 and candidate.orientation == "reverse")
            else "-"
        ),
        "source": source,
        "length_nt": len(candidate.sequence),
        "tm_c": candidate.tm_c,
        "gc_percent": candidate.gc_percent,
        "hairpin_tm_c": candidate.hairpin_tm_c,
        "homodimer_tm_c": candidate.homodimer_tm_c,
        "binding_start0": relative_start0,
        "binding_end0": relative_end0,
    }
    if region_interval0 is not None:
        result["genomic_interval_1based"] = _genomic_interval(
            region_start0=region_interval0[0],
            region_end0=region_interval0[1],
            offset_start0=relative_start0,
            offset_end0=relative_end0,
            gene_strand=gene_strand,
        )
        result["outside_homology_arm"] = True
        result["present_in_donor_plasmid"] = False
    elif source == "payload":
        result["payload_interval_1based"] = f"{relative_start0 + 1}-{relative_end0}"
        result["distance_from_5prime_junction_nt"] = relative_start0
        result["distance_from_3prime_junction_nt"] = payload_length - relative_end0
    return result


def design_genotyping_primers(
    *,
    assembly: str,
    chromosome: str,
    gene_strand: int,
    external_five_sequence: str,
    external_five_interval0: tuple[int, int],
    uha: str,
    wt_uha: str,
    removed_gene_sequence: str,
    dha: str,
    wt_dha: str,
    external_three_sequence: str,
    external_three_interval0: tuple[int, int],
    payload: str,
    assembled_plasmid: str,
) -> dict[str, Any]:
    """Design WT and two junction assays with genomic primers outside both arms."""
    result: dict[str, Any] = {
        "status": "PASS",
        "ruleset": PRIMER_RULESET,
        "payload_junction_standoff_nt": PAYLOAD_JUNCTION_STANDOFF_NT,
        "chromosome": chromosome,
        "assembly": assembly,
        "payload_sha256": sha256(payload.encode("ascii")).hexdigest(),
        "specificity_checks": {
            "donor_plasmid": "PASS for every released external genomic primer",
            "reference_genome": "NOT RUN; confirm with Primer-BLAST before ordering",
        },
        "assays": {},
        "warnings": [
            "Primer uniqueness is checked against the assembled donor plasmid, not genome-wide. "
            "Confirm genomic specificity with Primer-BLAST before ordering."
        ],
    }
    if len(payload) < PAYLOAD_JUNCTION_STANDOFF_NT + PRIMER_MIN_SIZE:
        result["status"] = "WARNING"
        result["warnings"].append(
            "Payload is too short to place a payload primer at least 150 bp from a junction."
        )
        return result

    assay_specs: list[dict[str, Any]] = []

    five_template = external_five_sequence + uha + payload
    five_left = _candidate_primers(
        five_template,
        region_start0=max(0, len(external_five_sequence) - 180),
        region_end0=len(external_five_sequence),
        orientation="forward",
        preferred_boundary0=len(external_five_sequence),
    )
    five_payload_start0 = len(external_five_sequence) + len(uha)
    five_right = _candidate_primers(
        five_template,
        region_start0=five_payload_start0 + PAYLOAD_JUNCTION_STANDOFF_NT,
        region_end0=min(len(five_template), five_payload_start0 + 380),
        orientation="reverse",
        preferred_boundary0=five_payload_start0 + PAYLOAD_JUNCTION_STANDOFF_NT,
    )
    assay_specs.append(
        {
            "key": "five_prime_junction",
            "template": five_template,
            "left": five_left,
            "right": five_right,
            "genomic_sides": {"left"},
            "target": len(uha) + PAYLOAD_JUNCTION_STANDOFF_NT + 100,
            "left_source": ("genomic_5prime_external", 0, external_five_interval0),
            "right_source": ("payload", five_payload_start0, None),
            "priority_side": "right",
        }
    )

    three_template = payload + dha + external_three_sequence
    three_left = _candidate_primers(
        three_template,
        region_start0=max(0, len(payload) - 380),
        region_end0=len(payload) - PAYLOAD_JUNCTION_STANDOFF_NT,
        orientation="forward",
        preferred_boundary0=len(payload) - PAYLOAD_JUNCTION_STANDOFF_NT,
    )
    three_external_start0 = len(payload) + len(dha)
    three_right = _candidate_primers(
        three_template,
        region_start0=three_external_start0,
        region_end0=min(len(three_template), three_external_start0 + 180),
        orientation="reverse",
        preferred_boundary0=three_external_start0,
    )
    assay_specs.append(
        {
            "key": "three_prime_junction",
            "template": three_template,
            "left": three_left,
            "right": three_right,
            "genomic_sides": {"right"},
            "target": len(dha) + PAYLOAD_JUNCTION_STANDOFF_NT + 100,
            "left_source": ("payload", 0, None),
            "right_source": (
                "genomic_3prime_external",
                three_external_start0,
                external_three_interval0,
            ),
            "priority_side": "left",
        }
    )

    wt_template = (
        external_five_sequence
        + wt_uha
        + removed_gene_sequence
        + wt_dha
        + external_three_sequence
    )
    wt_left = _candidate_primers(
        wt_template,
        region_start0=max(0, len(external_five_sequence) - 180),
        region_end0=len(external_five_sequence),
        orientation="forward",
        preferred_boundary0=len(external_five_sequence),
    )
    wt_external_start0 = (
        len(external_five_sequence)
        + len(wt_uha)
        + len(removed_gene_sequence)
        + len(wt_dha)
    )
    wt_right = _candidate_primers(
        wt_template,
        region_start0=wt_external_start0,
        region_end0=min(len(wt_template), wt_external_start0 + 180),
        orientation="reverse",
        preferred_boundary0=wt_external_start0,
    )
    assay_specs.append(
        {
            "key": "wild_type_locus",
            "template": wt_template,
            "left": wt_left,
            "right": wt_right,
            "genomic_sides": {"left", "right"},
            "target": len(wt_uha) + len(removed_gene_sequence) + len(wt_dha) + 80,
            "left_source": ("genomic_5prime_external", 0, external_five_interval0),
            "right_source": (
                "genomic_3prime_external",
                wt_external_start0,
                external_three_interval0,
            ),
            "priority_side": None,
        }
    )

    for spec in assay_specs:
        selected = _select_pair(
            assay=spec["key"],
            template=spec["template"],
            left_candidates=spec["left"],
            right_candidates=spec["right"],
            assembled_plasmid=assembled_plasmid,
            genomic_sides=spec["genomic_sides"],
            target_product_size=spec["target"],
            priority_side=spec["priority_side"],
        )
        if selected is None:
            result["status"] = "WARNING"
            result["assays"][spec["key"]] = {
                "assay": spec["key"],
                "status": "NO PRIMER PAIR",
                "reason": "No pair passed the configured sequence, Tm, structure, geometry, and donor-plasmid checks.",
            }
            continue
        left, right, assay = selected
        left_source, left_offset, left_interval = spec["left_source"]
        right_source, right_offset, right_interval = spec["right_source"]
        assay["forward_primer"] = _primer_dict(
            left,
            name=f"{spec['key']}_forward",
            source=left_source,
            source_offset0=left_offset,
            region_interval0=left_interval,
            gene_strand=gene_strand,
            payload_length=len(payload),
        )
        assay["reverse_primer"] = _primer_dict(
            right,
            name=f"{spec['key']}_reverse",
            source=right_source,
            source_offset0=right_offset,
            region_interval0=right_interval,
            gene_strand=gene_strand,
            payload_length=len(payload),
        )
        for primer in (assay["forward_primer"], assay["reverse_primer"]):
            if primer.get("source") == "payload":
                primer["reusable_payload_primer"] = True
                primer["reusable_payload_key"] = result["payload_sha256"]
        if spec["key"] == "wild_type_locus":
            edited_product_size = (
                int(assay["product_size_bp"])
                - len(removed_gene_sequence)
                + len(payload)
            )
            assay["allele_behavior"] = {
                "wild_type": "Expected to amplify",
                "edited": "Expected to amplify as a larger product if the PCR extension time supports it",
                "donor_plasmid": "Both external genomic primers are absent from the donor plasmid",
            }
            assay["expected_wild_type_product_size_bp"] = assay["product_size_bp"]
            assay["expected_edited_product_size_bp"] = edited_product_size
            assay["size_difference_bp"] = edited_product_size - int(assay["product_size_bp"])
            assay["expected_wild_type_amplicon_sequence_5to3"] = assay[
                "amplicon_sequence_5to3"
            ]
            edited_template = (
                external_five_sequence
                + uha
                + payload
                + dha
                + external_three_sequence
            )
            edited_right_end0 = right.bind_end0 + assay["size_difference_bp"]
            assay["expected_edited_amplicon_sequence_5to3"] = edited_template[
                left.bind_start0:edited_right_end0
            ]
        else:
            assay["allele_behavior"] = {
                "wild_type": "Payload primer is absent; no amplification expected",
                "edited": "Expected to amplify the verified insertion junction",
                "donor_plasmid": "External genomic primer is absent from the donor plasmid",
            }
        result["assays"][spec["key"]] = assay
    return result
