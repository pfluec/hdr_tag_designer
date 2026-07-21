from __future__ import annotations

import csv
import io

from .sequence import reverse_complement


def guide_oligos(guide_name: str, spacer: str) -> list[dict[str, str]]:
    """Return the two cloning oligos for one selected 20-nt guide spacer."""
    normalized_name = guide_name.strip()
    normalized_spacer = spacer.strip().upper()
    if not normalized_name:
        raise ValueError("Guide name must not be empty")
    if len(normalized_spacer) != 20 or set(normalized_spacer) - set("ACGT"):
        raise ValueError("Guide spacer must be exactly 20 A/C/G/T bases")
    return [
        {
            "Guide Name": normalized_name,
            "Sequence Type": f"{normalized_name}-fwd",
            "Sequence": "CACCG" + normalized_spacer,
        },
        {
            "Guide Name": normalized_name,
            "Sequence Type": f"{normalized_name}-rev",
            "Sequence": "AAAC" + reverse_complement(normalized_spacer) + "C",
        },
    ]


def guide_oligos_csv(guide_name: str, spacer: str) -> str:
    """Return the selected guide's two oligos in the supplied ordering CSV shape."""
    output = io.StringIO()
    writer = csv.DictWriter(
        output,
        fieldnames=["Guide Name", "Sequence Type", "Sequence"],
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerows(guide_oligos(guide_name, spacer))
    return output.getvalue()
