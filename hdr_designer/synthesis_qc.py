from __future__ import annotations

from typing import Any


TWIST_MAX_HOMOPOLYMER_NT = 14
TWIST_ORDERING_RULESET = "HDR Tag Designer Twist preflight v2 (2026-07-21)"
MIN_RETAINED_HOMOLOGY_ARM_NT = 100


def homopolymer_findings(
    sequence: str,
    *,
    max_allowed_nt: int = TWIST_MAX_HOMOPOLYMER_NT,
) -> list[dict[str, Any]]:
    """Return every homopolymer run longer than ``max_allowed_nt``."""
    normalized = sequence.upper()
    findings: list[dict[str, Any]] = []
    run_start0 = 0
    for position0 in range(1, len(normalized) + 1):
        if position0 < len(normalized) and normalized[position0] == normalized[run_start0]:
            continue
        run_length = position0 - run_start0
        if run_length > max_allowed_nt:
            findings.append(
                {
                    "base": normalized[run_start0],
                    "start0": run_start0,
                    "end0": position0,
                    "interval_1based": f"{run_start0 + 1}-{position0}",
                    "length_nt": run_length,
                }
            )
        run_start0 = position0
    return findings
