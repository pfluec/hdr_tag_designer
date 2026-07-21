from __future__ import annotations

import copy
import csv
import hashlib
import io
import unittest
import zipfile

from Bio import SeqIO

from hdr_designer.design import (
    adjust_homology_arm_boundary_for_synthesis,
    design_tubb5_fixture,
)
from hdr_designer.exports import genotyping_primers_csv
from hdr_designer.generate_oligos_from_guides import guide_oligos, guide_oligos_csv
from hdr_designer.models import HomologyArm
from hdr_designer.ordering import (
    OrderingError,
    design_identity,
    homopolymer_findings,
    ordering_package_members,
    ordering_package_zip,
    twist_ordering_qc,
    twist_sequences_csv,
)
from hdr_designer.sequence import gc_percent, reverse_complement


class OrderingExportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tubb5 = design_tubb5_fixture()
        cls.ordering_ready = cls.tubb5

    @staticmethod
    def _arm(sequence: str, strand: int, name: str) -> HomologyArm:
        return HomologyArm(
            name=name,
            length=len(sequence),
            chromosome="1",
            genomic_start0=1_000,
            genomic_end0=1_000 + len(sequence),
            gene_oriented_sequence=sequence,
            chromosome_forward_sequence=(
                sequence if strand == 1 else reverse_complement(sequence)
            ),
            gc_percent=gc_percent(sequence),
            requested_length=len(sequence),
        )

    def test_selected_guide_oligos_follow_supplied_function(self) -> None:
        rows = guide_oligos("Tubb5_selected", "GAGGCAGAAGAGGAGGCCTA")
        self.assertEqual(
            rows,
            [
                {
                    "Guide Name": "Tubb5_selected",
                    "Sequence Type": "Tubb5_selected-fwd",
                    "Sequence": "CACCGGAGGCAGAAGAGGAGGCCTA",
                },
                {
                    "Guide Name": "Tubb5_selected",
                    "Sequence Type": "Tubb5_selected-rev",
                    "Sequence": "AAACTAGGCCTCCTCTTCTGCCTCC",
                },
            ],
        )
        parsed = list(csv.DictReader(io.StringIO(guide_oligos_csv("Tubb5_selected", rows[0]["Sequence"][5:]))))
        self.assertEqual(parsed, rows)
        with self.assertRaisesRegex(ValueError, "exactly 20"):
            guide_oligos("bad", "ACGT")

    def test_genotyping_csv_contains_four_unique_shared_primers(self) -> None:
        rows = list(csv.DictReader(io.StringIO(genotyping_primers_csv(self.tubb5))))
        self.assertEqual(
            [row["primer_name"] for row in rows],
            [
                "Tubb5_wt_5_fwd",
                "Tubb5_wt_3_rev",
                "Tubb5_mut_5_rev",
                "Tubb5_mut_3_fwd",
            ],
        )
        self.assertEqual(len({row["sequence_5to3"] for row in rows}), 4)
        assays = self.tubb5.genotyping_primers["assays"]
        self.assertEqual(
            assays["wild_type_locus"]["forward_primer"]["sequence_5to3"],
            assays["five_prime_junction"]["forward_primer"]["sequence_5to3"],
        )
        self.assertEqual(
            assays["wild_type_locus"]["reverse_primer"]["sequence_5to3"],
            assays["three_prime_junction"]["reverse_primer"]["sequence_5to3"],
        )

    def test_homopolymer_limit_is_strictly_greater_than_fourteen(self) -> None:
        self.assertEqual(homopolymer_findings("C" + "A" * 14 + "G"), [])
        self.assertEqual(
            homopolymer_findings("C" + "A" * 15 + "G"),
            [
                {
                    "base": "A",
                    "start0": 1,
                    "end0": 16,
                    "interval_1based": "2-16",
                    "length_nt": 15,
                }
            ],
        )

    def test_uha_boundary_keeps_insertion_proximal_suffix_on_both_strands(self) -> None:
        sequence = "CG" * 100 + "A" * 20 + "TG" * 190
        for strand in (1, -1):
            with self.subTest(strand=strand):
                original = self._arm(sequence, strand, "5-prime homology arm")
                adjusted, note = adjust_homology_arm_boundary_for_synthesis(
                    original,
                    arm_role="five_prime",
                    gene_strand=strand,
                )
                self.assertEqual(adjusted.length, 390)
                self.assertTrue(adjusted.gene_oriented_sequence.startswith("A" * 10))
                self.assertEqual(homopolymer_findings(adjusted.gene_oriented_sequence), [])
                self.assertIn("guide/SapI-facing arm edge", str(note))
                if strand == 1:
                    self.assertEqual(
                        (adjusted.genomic_start0, adjusted.genomic_end0),
                        (1_210, 1_600),
                    )
                else:
                    self.assertEqual(
                        (adjusted.genomic_start0, adjusted.genomic_end0),
                        (1_000, 1_390),
                    )

    def test_dha_boundary_keeps_insertion_proximal_prefix_on_both_strands(self) -> None:
        sequence = "CG" * 150 + "T" * 20 + "AG" * 140
        for strand in (1, -1):
            with self.subTest(strand=strand):
                original = self._arm(sequence, strand, "3-prime homology arm")
                adjusted, _ = adjust_homology_arm_boundary_for_synthesis(
                    original,
                    arm_role="three_prime",
                    gene_strand=strand,
                )
                self.assertEqual(adjusted.length, 310)
                self.assertTrue(adjusted.gene_oriented_sequence.endswith("T" * 10))
                self.assertEqual(homopolymer_findings(adjusted.gene_oriented_sequence), [])
                if strand == 1:
                    self.assertEqual(
                        (adjusted.genomic_start0, adjusted.genomic_end0),
                        (1_000, 1_310),
                    )
                else:
                    self.assertEqual(
                        (adjusted.genomic_start0, adjusted.genomic_end0),
                        (1_290, 1_600),
                    )

    def test_multiple_runs_select_boundary_that_removes_all_long_runs(self) -> None:
        sequence = (
            "AG" * 60
            + "A" * 20
            + "TG" * 80
            + "C" * 18
            + "AT" * 141
        )
        uha, _ = adjust_homology_arm_boundary_for_synthesis(
            self._arm(sequence, 1, "5-prime homology arm"),
            arm_role="five_prime",
            gene_strand=1,
        )
        dha, _ = adjust_homology_arm_boundary_for_synthesis(
            self._arm(sequence, 1, "3-prime homology arm"),
            arm_role="three_prime",
            gene_strand=1,
        )
        self.assertTrue(uha.gene_oriented_sequence.startswith("C" * 9))
        self.assertTrue(dha.gene_oriented_sequence.endswith("A" * 10))
        self.assertEqual(homopolymer_findings(uha.gene_oriented_sequence), [])
        self.assertEqual(homopolymer_findings(dha.gene_oriented_sequence), [])

    def test_tubb5_moves_dha_boundary_into_homopolymer(self) -> None:
        qc = twist_ordering_qc(self.tubb5)
        self.assertEqual(qc["status"], "PASS")
        self.assertEqual(qc["findings"], [])
        adjustment = self.tubb5.three_prime_arm.boundary_adjustment
        self.assertEqual(
            adjustment,
            {
                "status": "ADJUSTED",
                "reason": "Twist homopolymer ordering limit",
                "rule_max_nt": 14,
                "original_length_nt": 600,
                "final_length_nt": 268,
                "trimmed_bases_nt": 332,
                "original_genomic_interval_1based": "17:36145274-36145873",
                "final_genomic_interval_1based": "17:36145606-36145873",
                "boundary_side": "gene-oriented 3-prime distal boundary",
                "homopolymer_base": "T",
                "original_run_interval_1based": "256-281",
                "original_run_length_nt": 26,
                "retained_boundary_run_length_nt": 13,
            },
        )
        self.assertTrue(self.tubb5.three_prime_arm.final_gene_oriented_sequence.endswith("T" * 13))
        self.assertEqual(qc["boundary_adjustments"][0]["arm"], "3-prime homology arm")
        self.assertGreater(len(ordering_package_zip(self.tubb5)), 0)

    def test_unresolved_final_arm_homopolymer_remains_a_package_guard(self) -> None:
        blocked = copy.deepcopy(self.tubb5)
        sequence = blocked.three_prime_arm.final_gene_oriented_sequence
        blocked.three_prime_arm.corrected_gene_oriented_sequence = "A" * 15 + sequence[15:]
        self.assertEqual(twist_ordering_qc(blocked)["status"], "ERROR")
        with self.assertRaisesRegex(OrderingError, "homopolymer validation failed"):
            ordering_package_zip(blocked)

    def test_twist_csv_contains_final_synthesis_fragments(self) -> None:
        rows = list(csv.DictReader(io.StringIO(twist_sequences_csv(self.ordering_ready))))
        self.assertEqual(len(rows), 2)
        expected_sequences = [
            self.ordering_ready.cloning_fragments["uha_synthesis_fragment_5to3"],
            self.ordering_ready.cloning_fragments["dha_synthesis_fragment_5to3"],
        ]
        self.assertEqual([row["Sequence"] for row in rows], expected_sequences)
        self.assertEqual([row["Sequence Type"] for row in rows], [
            "final_uha_synthesis_fragment",
            "final_dha_synthesis_fragment",
        ])
        for row, sequence in zip(rows, expected_sequences):
            self.assertEqual(int(row["Length"]), len(sequence))
            self.assertEqual(row["SHA-256"], hashlib.sha256(sequence.encode("ascii")).hexdigest())
            self.assertEqual(row["Internal QC Status"], "PASS")
            self.assertEqual(row["Internal QC Ruleset"], "HDR Tag Designer Twist preflight v2 (2026-07-21)")
            self.assertEqual(row["Twist Portal Screening"], "REQUIRED")
        self.assertEqual(rows[0]["Requested Arm Length"], "600")
        self.assertEqual(rows[0]["Final Arm Length"], "600")
        self.assertEqual(rows[1]["Requested Arm Length"], "600")
        self.assertEqual(rows[1]["Final Arm Length"], "268")
        self.assertIn("shortened from 600 to 268 bp", rows[1]["Boundary Adjustment"])

    def test_zip_is_deterministic_and_contains_only_supported_outputs(self) -> None:
        first = ordering_package_zip(self.ordering_ready)
        second = ordering_package_zip(self.ordering_ready)
        self.assertEqual(first, second)
        expected_names = [name for name, _ in ordering_package_members(self.ordering_ready)]
        stem = f"Tubb5_c_terminal_{design_identity(self.ordering_ready)}"
        self.assertEqual(
            expected_names,
            [
                f"{stem}_twist_sequences.csv",
                f"{stem}_guide_oligos.csv",
                f"{stem}_genotyping_primers.csv",
                f"{stem}_assembled_plasmid.gb",
                f"{stem}_wild_type_locus.gb",
                f"{stem}_edited_locus.gb",
            ],
        )
        with zipfile.ZipFile(io.BytesIO(first)) as archive:
            self.assertEqual(archive.namelist(), expected_names)
            self.assertTrue(
                all(item.date_time == (1980, 1, 1, 0, 0, 0) for item in archive.infolist())
            )
            guide_rows = list(
                csv.DictReader(io.StringIO(archive.read(expected_names[1]).decode("utf-8")))
            )
            self.assertEqual(len(guide_rows), 2)
            primer_rows = list(
                csv.DictReader(io.StringIO(archive.read(expected_names[2]).decode("utf-8")))
            )
            self.assertEqual(len(primer_rows), 4)
            for name in expected_names[3:]:
                record = SeqIO.read(io.StringIO(archive.read(name).decode("utf-8")), "genbank")
                self.assertGreater(len(record.seq), 0)
                self.assertEqual(record.annotations["date"], "01-JAN-1980")


if __name__ == "__main__":
    unittest.main()
