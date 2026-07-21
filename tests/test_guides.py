from __future__ import annotations

import unittest

from hdr_designer.guides import enumerate_spcas9_guides


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


if __name__ == "__main__":
    unittest.main()
