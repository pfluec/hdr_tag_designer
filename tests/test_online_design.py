from __future__ import annotations

import unittest

from hdr_designer.design import design_online
from hdr_designer.ensembl import SPECIES
from hdr_designer.models import Exon, TranscriptRecord


class SyntheticEnsemblClient:
    """Deterministic Ensembl stand-in for exercising live-design code paths."""

    def __init__(self, species_key: str) -> None:
        self.species = SPECIES[species_key]
        chromosome = list("A" * 4000)
        exon_start0 = 1000
        exon_end0 = 2500
        cdna = list("A" * (exon_end0 - exon_start0))
        cds_start0 = 31
        c_terminal_boundary0 = 2000
        stop_cdna0 = c_terminal_boundary0 - exon_start0

        coding_before_stop = "ATG" + "GCT" * (
            (stop_cdna0 - cds_start0 - 3) // 3
        )
        self.assert_equal_length(coding_before_stop, stop_cdna0 - cds_start0)
        cdna[cds_start0:stop_cdna0] = coding_before_stop
        cdna[stop_cdna0:stop_cdna0 + 3] = "TAA"

        # One safe NGG candidate near each terminus. Both targets are split by
        # the intended edit and contain no extra nearby NGG or CCN motifs.
        targets = [
            (
                exon_start0 + cds_start0 + 3 - 10,
                "ACACACAATGACACACACACAGG",
            ),
            (
                c_terminal_boundary0 - 10,
                "ACACACACACTAATACACACAGG",
            ),
        ]
        for start0, target in targets:
            self.assert_equal_length(target, 23)
            chromosome[start0:start0 + 23] = target
            local_start0 = start0 - exon_start0
            cdna[local_start0:local_start0 + 23] = target

        cdna_text = "".join(cdna)
        cds = cdna_text[cds_start0:stop_cdna0 + 3]
        if not cds.startswith("ATG") or not cds.endswith("TAA"):
            raise AssertionError("Synthetic CDS boundary construction failed")
        if len(cds[:-3]) % 3:
            raise AssertionError("Synthetic CDS is out of frame")

        self.record = TranscriptRecord(
            species=self.species,
            gene_symbol="MockTagGene",
            gene_id="MOCKG0001",
            transcript_id="MOCKT0001",
            transcript_version="1",
            chromosome="1",
            strand=1,
            cdna=cdna_text,
            cds=cds,
            exons=[
                Exon(
                    start1=exon_start0 + 1,
                    end1=exon_end0,
                    strand=1,
                    stable_id="MOCKE1",
                )
            ],
            source="synthetic test fixture",
            source_release="test",
        )
        self.chromosome = "".join(chromosome)

    @staticmethod
    def assert_equal_length(value: str, expected: int) -> None:
        if len(value) != expected:
            raise AssertionError(f"Expected {expected} nt, got {len(value)}")

    def transcript_record(
        self,
        species,
        symbol_or_gene_id: str,
        transcript_id: str | None = None,
    ) -> TranscriptRecord:
        if species != self.species:
            raise AssertionError("Wrong species supplied to synthetic client")
        return self.record

    def region_sequence(
        self,
        species,
        chromosome: str,
        start0: int,
        end0: int,
    ) -> str:
        if species != self.species or chromosome != "1":
            raise AssertionError("Wrong region supplied to synthetic client")
        return self.chromosome[start0:end0]


class OnlineDesignPathTest(unittest.TestCase):
    def test_mouse_and_human_c_terminal_paths(self) -> None:
        for species_key, assembly in (("mouse", "GRCm39"), ("human", "GRCh38")):
            with self.subTest(species=species_key):
                result = design_online(
                    species_key=species_key,
                    gene="MockTagGene",
                    client=SyntheticEnsemblClient(species_key),
                )
                self.assertTrue(result.sequence_complete)
                self.assertEqual(result.assembly, assembly)
                self.assertEqual(result.top_guide.distance_to_insertion, 7)
                self.assertFalse(result.top_guide.blocking_mutation_required)
                self.assertEqual(
                    result.cloning_fragments["assembled_plasmid_length_nt"],
                    3950,
                )
                self.assertEqual(
                    result.cloning_fragments["assembled_plasmid_sapi_site_count"],
                    0,
                )

    def test_n_terminal_path_is_explicit_preview(self) -> None:
        result = design_online(
            species_key="human",
            gene="MockTagGene",
            terminus="N-terminal",
            client=SyntheticEnsemblClient("human"),
        )
        self.assertFalse(result.sequence_complete)
        self.assertEqual(result.status, "REVIEW REQUIRED")
        self.assertIn("uha_synthesis_fragment_preview_5to3", result.cloning_fragments)
        self.assertNotIn("assembled_plasmid_5to3", result.cloning_fragments)
        self.assertTrue(
            any("N-terminal output is a locus-design preview" in warning for warning in result.warnings)
        )


if __name__ == "__main__":
    unittest.main()
