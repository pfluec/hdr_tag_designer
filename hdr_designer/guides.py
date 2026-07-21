from __future__ import annotations

from .models import GuideCandidate
from .sequence import gc_percent, reverse_complement

GUIDE_SAFETY_CUTOFF_NT = 14


def _interval_overlap(start_a: int, end_a: int, start_b: int, end_b: int) -> bool:
    return start_a < end_b and start_b < end_a


def _insertion_splits_interval(boundary: int, start: int, end: int) -> bool:
    return start < boundary < end


def _retained_segment_length(
    target_start0: int,
    target_end0: int,
    insertion_boundary0: int,
    removed_start0: int,
    removed_end0: int,
) -> int:
    """Longest contiguous part of the original 23-nt target after the intended edit."""
    segments: list[tuple[int, int]] = [(target_start0, target_end0)]

    if removed_start0 < removed_end0:
        after_deletion: list[tuple[int, int]] = []
        for start0, end0 in segments:
            if not _interval_overlap(start0, end0, removed_start0, removed_end0):
                after_deletion.append((start0, end0))
                continue
            if start0 < removed_start0:
                after_deletion.append((start0, min(end0, removed_start0)))
            if removed_end0 < end0:
                after_deletion.append((max(start0, removed_end0), end0))
        segments = after_deletion

    after_insertion: list[tuple[int, int]] = []
    for start0, end0 in segments:
        if start0 < insertion_boundary0 < end0:
            after_insertion.extend(
                [(start0, insertion_boundary0), (insertion_boundary0, end0)]
            )
        else:
            after_insertion.append((start0, end0))

    return max((end0 - start0 for start0, end0 in after_insertion), default=0)


def _activity_heuristic(spacer: str) -> str:
    gc = gc_percent(spacer)
    if "TTTT" in spacer:
        return "unfavorable: contains TTTT"
    if 40.0 <= gc <= 60.0:
        return "favorable basic spacer properties"
    if 30.0 <= gc <= 70.0:
        return "acceptable basic spacer properties"
    return "review: extreme GC content"



def longest_retained_segment_after_point_mutations(
    *,
    target_start0: int,
    target_end0: int,
    insertion_boundary0: int,
    removed_start0: int,
    removed_end0: int,
    mutated_genomic_positions0: list[int],
) -> int:
    """Longest original target segment after edit boundaries and point mismatches."""
    segments: list[tuple[int, int]] = [(target_start0, target_end0)]
    if removed_start0 < removed_end0:
        trimmed: list[tuple[int, int]] = []
        for start0, end0 in segments:
            if not _interval_overlap(start0, end0, removed_start0, removed_end0):
                trimmed.append((start0, end0))
            else:
                if start0 < removed_start0:
                    trimmed.append((start0, min(end0, removed_start0)))
                if removed_end0 < end0:
                    trimmed.append((max(start0, removed_end0), end0))
        segments = trimmed
    split: list[tuple[int, int]] = []
    for start0, end0 in segments:
        if start0 < insertion_boundary0 < end0:
            split.extend([(start0, insertion_boundary0), (insertion_boundary0, end0)])
        else:
            split.append((start0, end0))
    segments = split
    for position0 in sorted(set(mutated_genomic_positions0)):
        split = []
        for start0, end0 in segments:
            if start0 <= position0 < end0:
                if start0 < position0:
                    split.append((start0, position0))
                if position0 + 1 < end0:
                    split.append((position0 + 1, end0))
            else:
                split.append((start0, end0))
        segments = split
    return max((end0 - start0 for start0, end0 in segments), default=0)

def enumerate_spcas9_guides(
    chromosome_forward_sequence: str,
    region_start0: int,
    insertion_boundary0: int,
    search_radius: int,
    removed_start0: int,
    removed_end0: int,
) -> list[GuideCandidate]:
    """Enumerate SpCas9-NGG sites with a nominal nick inside the requested window.

    Coordinates are 0-based half-open genomic coordinates. ``spacer`` and ``pam``
    are reported 5' to 3' on the PAM-containing strand. The rank follows the Bollen
    protocol priorities available in this prototype: nick proximity first, whether
    the intended edit disrupts the target second, then basic spacer heuristics.
    """
    sequence = chromosome_forward_sequence.upper()
    region_end0 = region_start0 + len(sequence)
    desired_start0 = insertion_boundary0 - search_radius
    desired_end0 = insertion_boundary0 + search_radius
    candidates: list[GuideCandidate] = []

    for local_start in range(0, max(0, len(sequence) - 22)):
        genomic_start0 = region_start0 + local_start
        target_forward = sequence[local_start:local_start + 23]
        if len(target_forward) != 23:
            continue

        # PAM on chromosome-forward strand: protospacer-NGG.
        if target_forward[21:23] == "GG":
            nick_boundary0 = genomic_start0 + 17
            if desired_start0 <= nick_boundary0 <= desired_end0:
                candidates.append(
                    _make_candidate(
                        spacer=target_forward[:20],
                        pam=target_forward[20:23],
                        chromosome_strand="+",
                        target_start0=genomic_start0,
                        target_end0=genomic_start0 + 23,
                        pam_start0=genomic_start0 + 20,
                        pam_end0=genomic_start0 + 23,
                        nick_boundary0=nick_boundary0,
                        insertion_boundary0=insertion_boundary0,
                        removed_start0=removed_start0,
                        removed_end0=removed_end0,
                    )
                )

        # PAM on chromosome-reverse strand appears as CCN on forward sequence.
        if target_forward[:2] == "CC":
            nick_boundary0 = genomic_start0 + 6
            if desired_start0 <= nick_boundary0 <= desired_end0:
                candidates.append(
                    _make_candidate(
                        spacer=reverse_complement(target_forward[3:23]),
                        pam=reverse_complement(target_forward[:3]),
                        chromosome_strand="-",
                        target_start0=genomic_start0,
                        target_end0=genomic_start0 + 23,
                        pam_start0=genomic_start0,
                        pam_end0=genomic_start0 + 3,
                        nick_boundary0=nick_boundary0,
                        insertion_boundary0=insertion_boundary0,
                        removed_start0=removed_start0,
                        removed_end0=removed_end0,
                    )
                )

    candidates = [
        candidate
        for candidate in candidates
        if region_start0 <= candidate.target_start0
        and candidate.target_end0 <= region_end0
    ]
    candidates.sort(
        key=lambda candidate: (
            candidate.distance_to_insertion,
            not candidate.target_destroyed,
            candidate.poly_t,
            abs(candidate.gc_percent - 50.0),
            candidate.spacer,
        )
    )
    for rank, candidate in enumerate(candidates, start=1):
        candidate.rank = rank
    return candidates


def _make_candidate(
    *,
    spacer: str,
    pam: str,
    chromosome_strand: str,
    target_start0: int,
    target_end0: int,
    pam_start0: int,
    pam_end0: int,
    nick_boundary0: int,
    insertion_boundary0: int,
    removed_start0: int,
    removed_end0: int,
) -> GuideCandidate:
    insertion_disrupts_target = _insertion_splits_interval(
        insertion_boundary0, target_start0, target_end0
    )
    deletion_disrupts_target = _interval_overlap(
        target_start0, target_end0, removed_start0, removed_end0
    )
    target_destroyed = insertion_disrupts_target or deletion_disrupts_target
    pam_destroyed = _insertion_splits_interval(
        insertion_boundary0, pam_start0, pam_end0
    ) or _interval_overlap(pam_start0, pam_end0, removed_start0, removed_end0)
    retained = _retained_segment_length(
        target_start0,
        target_end0,
        insertion_boundary0,
        removed_start0,
        removed_end0,
    )
    gc = gc_percent(spacer)
    poly_t = "TTTT" in spacer
    # Bollen's 14-nt safeguard applies to a surviving PAM-proximal target. If
    # the intended edit removes part of the PAM itself, the edited allele no
    # longer contains a functional SpCas9-NGG target and no extra blocking
    # mutation is required.
    blocking_required = (not pam_destroyed) and retained > GUIDE_SAFETY_CUTOFF_NT

    reasons = [f"nominal nick {abs(nick_boundary0 - insertion_boundary0)} bp from insertion"]
    if pam_destroyed:
        reasons.append("PAM overlaps the intended edit and is destroyed")
    elif target_destroyed:
        reasons.append("protospacer is split by the intended edit")
    else:
        reasons.append("target is not disrupted by the intended edit")
    reasons.append(f"longest retained target segment {retained} nt")
    if pam_destroyed:
        reasons.append("no additional blocking mutation required because the PAM is lost")
    elif blocking_required:
        reasons.append(f">{GUIDE_SAFETY_CUTOFF_NT} nt with intact PAM: blocking mutation required")
    else:
        reasons.append(f"within the {GUIDE_SAFETY_CUTOFF_NT}-nt protocol cutoff")
    reasons.append(f"GC {gc:.1f}%")
    if poly_t:
        reasons.append("contains TTTT")

    if pam_destroyed:
        blocking_note = (
            "No extra guide-blocking mutation is required because the intended edit "
            "removes part of the NGG PAM."
        )
    elif blocking_required:
        blocking_note = (
            "A guide-blocking change is required before order-ready donor sequences can "
            f"be released: the PAM remains intact and {retained} nt of the original "
            f"target remain contiguous (protocol cutoff <= {GUIDE_SAFETY_CUTOFF_NT} nt)."
        )
    else:
        blocking_note = (
            "No extra guide-blocking mutation is required under the protocol cutoff."
        )

    return GuideCandidate(
        rank=0,
        spacer=spacer,
        pam=pam,
        chromosome_strand=chromosome_strand,
        target_start0=target_start0,
        target_end0=target_end0,
        pam_start0=pam_start0,
        pam_end0=pam_end0,
        nick_boundary0=nick_boundary0,
        distance_to_insertion=abs(nick_boundary0 - insertion_boundary0),
        gc_percent=gc,
        poly_t=poly_t,
        insertion_disrupts_target=insertion_disrupts_target,
        deletion_disrupts_target=deletion_disrupts_target,
        target_destroyed=target_destroyed,
        pam_destroyed=pam_destroyed,
        longest_retained_segment=retained,
        final_longest_retained_segment=retained,
        blocking_mutation_required=blocking_required,
        activity_heuristic=_activity_heuristic(spacer),
        rationale="; ".join(reasons),
        blocking_mutation_note=blocking_note,
        final_pam=pam,
        final_pam_destroyed=pam_destroyed,
    )
