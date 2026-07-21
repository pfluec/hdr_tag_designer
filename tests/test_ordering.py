from __future__ import annotations

import copy
import csv
import hashlib
import io
import unittest
import zipfile

from Bio import SeqIO

from hdr_designer.design import design_tubb5_fixture
from hdr_designer.exports import genotyping_primers_csv
from hdr_designer.generate_oligos_from_guides import guide_oligos, guide_oligos_csv
from hdr_designer.ordering import (
    OrderingError,
    design_identity,
    homopolymer_findings,
    ordering_package_members,
    ordering_package_zip,
    twist_ordering_qc,
    twist_sequences_csv,
)


class OrderingExportTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.tubb5 = design_tubb5_fixture()
        cls.ordering_ready = copy.deepcopy(cls.tubb5)

        arm = cls.ordering_ready.three_prime_arm
        old_arm = arm.final_gene_oriented_sequence
        findings = homopolymer_findings(old_arm)
        assert len(findings) == 1
        position0 = findings[0]["start0"] + findings[0]["length_nt"] // 2
        replacement = "A" if old_arm[position0] != "A" else "C"
        new_arm = old_arm[:position0] + replacement + old_arm[position0 + 1 :]
        arm.corrected_gene_oriented_sequence = new_arm

        for key in (
            "dha_synthesis_fragment_5to3",
            "assembled_donor_insert_5to3",
            "assembled_plasmid_5to3",
        ):
            value = str(cls.ordering_ready.cloning_fragments[key])
            cls.ordering_ready.cloning_fragments[key] = value.replace(old_arm, new_arm, 1)
        edited = cls.ordering_ready.locus_contexts["edited"]
        edited["sequence_5to3"] = str(edited["sequence_5to3"]).replace(old_arm, new_arm, 1)

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

    def test_tubb5_is_blocked_by_exact_final_arm_homopolymer(self) -> None:
        qc = twist_ordering_qc(self.tubb5)
        self.assertEqual(qc["status"], "ERROR")
        self.assertEqual(
            qc["findings"],
            [
                {
                    "arm": "3-prime homology arm",
                    "base": "T",
                    "start0": 255,
                    "end0": 281,
                    "interval_1based": "256-281",
                    "length_nt": 26,
                }
            ],
        )
        with self.assertRaisesRegex(OrderingError, "T x 26"):
            ordering_package_zip(self.tubb5)

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
            self.assertEqual(row["Twist Portal Screening"], "REQUIRED")

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
