from __future__ import annotations

import unittest

from hdr_designer.guides import enumerate_spcas9_guides, recuttable_on_target_sites
from hdr_designer.sequence import reverse_complement


class GuidePriorityTest(unittest.TestCase):
    def test_nearer_guide_outranks_farther_pam_destroyed_guide(self) -> None:
        sequence = list("A" * 100)
        # Nick 6 bp from the insertion; the target is retained and needs blocking.
        sequence[27:50] = "ACACACACACACACACACACAGG"
        # Reverse-strand guide nick 7 bp away; its PAM overlaps the removed interval.
        sequence[51:74] = "CCAATATATATATATATATATAT"
        guides = enumerate_spcas9_guides(
            chromosome_forward_sequence="".join(sequence),
            region_start0=0,
            insertion_boundary0=50,
            search_radius=20,
            removed_start0=50,
            removed_end0=53,
        )
        self.assertGreaterEqual(len(guides), 2)
        self.assertEqual(guides[0].distance_to_insertion, 6)
        self.assertTrue(guides[0].blocking_mutation_required)
        farther_safe = next(
            guide
            for guide in guides
            if guide.distance_to_insertion == 7 and guide.pam_destroyed
        )
        self.assertGreater(farther_safe.rank, guides[0].rank)


class EditedLocusRecuttingTest(unittest.TestCase):
    spacer = "TCATGTCCTTACTGCTGTCA"

    def matches(self, sequence: str, *junctions: int):
        return recuttable_on_target_sites(
            sequence, self.spacer, junctions0=tuple(junctions)
        )

    def test_hck_insertion_separates_protospacer_and_pam(self) -> None:
        insertion = "ACGT" * 182 + "A"
        edited = self.spacer[:17] + insertion + "AGG"
        self.assertEqual(len(insertion), 729)
        self.assertEqual(self.matches(edited, 17, 746), [])

    def test_entire_spacer_and_pam_remain_contiguous(self) -> None:
        matches = self.matches(self.spacer + "AGG" + "A" * 20, 30)
        self.assertTrue(any(match["exact_original_site"] for match in matches))

    def test_destroyed_pam_is_not_a_candidate(self) -> None:
        self.assertEqual(self.matches(self.spacer + "AGA"), [])

    def test_target_recreated_across_insertion_junction(self) -> None:
        edited = "AAAA" + self.spacer[:10] + self.spacer[10:] + "TGG" + "AAAA"
        matches = self.matches(edited, 14)
        self.assertTrue(any(match["crosses_insertion_junction"] for match in matches))

    def test_reverse_complement_site_is_detected(self) -> None:
        edited = "AAAA" + reverse_complement(self.spacer + "CGG") + "AAAA"
        matches = self.matches(edited)
        self.assertTrue(any(match["strand"] == "-" for match in matches))


if __name__ == "__main__":
    unittest.main()
