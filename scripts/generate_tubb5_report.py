from __future__ import annotations

from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from hdr_designer.design import design_tubb5_fixture
from hdr_designer.exports import (
    arms_fasta,
    assembled_plasmid_genbank,
    design_json,
    design_report,
    guides_csv,
)


def main() -> None:
    output = ROOT / "outputs" / "tubb5_test"
    output.mkdir(parents=True, exist_ok=True)
    result = design_tubb5_fixture()
    (output / "Tubb5_C_terminal_design_report.txt").write_text(
        design_report(result), encoding="utf-8"
    )
    (output / "Tubb5_guides.csv").write_text(guides_csv(result), encoding="utf-8")
    (output / "Tubb5_homology_arms.fasta").write_text(
        arms_fasta(result), encoding="utf-8"
    )
    (output / "Tubb5_design.json").write_text(design_json(result), encoding="utf-8")
    (output / "Tubb5_C_terminal_mNeonGreen_assembled_plasmid.gb").write_text(
        assembled_plasmid_genbank(result), encoding="utf-8"
    )
    print(f"Wrote Tubb5 test outputs to {output}")


if __name__ == "__main__":
    main()
