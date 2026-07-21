from __future__ import annotations

from io import StringIO
import json
import unittest

from Bio import SeqIO

from hdr_designer.backbones import (
    BACKBONE_DNA,
    NTERM_BACKBONE_DNA,
    NTERM_DHA_PREFIX,
    NTERM_DHA_SUFFIX,
    NTERM_UHA_PREFIX,
    NTERM_UHA_SUFFIX,
    backbone_metadata_for,
    infer_custom_backbone_definition,
    n_terminal_backbone_analysis,
    NTERM_BACKBONE,
)
from hdr_designer.design import _automatic_noncoding_sapi_mutation, design_online
from hdr_designer.ensembl import SPECIES
from hdr_designer.exports import (
    assembled_plasmid_genbank,
    design_json,
    design_report,
    sapi_qc_rows,
)
from hdr_designer.models import Exon, HomologyArm, TranscriptRecord
from hdr_designer.sequence import reverse_complement, translate


class SyntheticEnsemblClient:
    """Deterministic Ensembl stand-in for exercising live-design code paths."""

    def __init__(
        self,
        species_key: str,
        *,
        guide_requires_blocking: bool = False,
        guide_requires_seed_blocking: bool = False,
        guide_has_no_synonymous_block: bool = False,
        coding_sapi_site: bool = False,
        multiple_coding_sapi_sites: bool = False,
        noncoding_sapi_site: bool = False,
        strand: int = 1,
    ) -> None:
        if strand not in {-1, 1}:
            raise ValueError("strand must be -1 or 1")
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
        if guide_requires_seed_blocking:
            c_terminal_target = (
                c_terminal_boundary0 - 23,
                "ACACACACACACACACACACTGG",
            )
        elif guide_has_no_synonymous_block:
            c_terminal_target = (
                c_terminal_boundary0 + 3,
                "ACACACACACACACACACACAGG",
            )
        elif guide_requires_blocking:
            c_terminal_target = (
                c_terminal_boundary0 - 23,
                "GCTGCTGCTGCTGCTGCTGCGGG",
            )
        else:
            c_terminal_target = (
                c_terminal_boundary0 - 10,
                "ACACACACACTAATACACACAGG",
            )
        targets = [
            (
                exon_start0 + cds_start0 + 3 - 10,
                "ACACACAATGACACACACACAGG",
            ),
            c_terminal_target,
        ]
        for start0, target in targets:
            self.assert_equal_length(target, 23)
            chromosome[start0:start0 + 23] = target
            local_start0 = start0 - exon_start0
            cdna[local_start0:local_start0 + 23] = target

        if coding_sapi_site:
            # GAG|CTC|TTC is coding and contains GCTCTTC across codons.
            sapi_cdna_start0 = stop_cdna0 - 300
            sapi_sequence = "GAGCTCTTC"
            cdna[sapi_cdna_start0:sapi_cdna_start0 + len(sapi_sequence)] = sapi_sequence
            sapi_genomic_start0 = exon_start0 + sapi_cdna_start0
            chromosome[sapi_genomic_start0:sapi_genomic_start0 + len(sapi_sequence)] = sapi_sequence

        if multiple_coding_sapi_sites:
            for offset in (300, 450):
                sapi_cdna_start0 = stop_cdna0 - offset
                sapi_sequence = "GAGCTCTTC"
                cdna[sapi_cdna_start0:sapi_cdna_start0 + len(sapi_sequence)] = sapi_sequence
                sapi_genomic_start0 = exon_start0 + sapi_cdna_start0
                chromosome[sapi_genomic_start0:sapi_genomic_start0 + len(sapi_sequence)] = sapi_sequence

        if noncoding_sapi_site:
            sapi_cdna_start0 = stop_cdna0 + 3 + 100
            sapi_sequence = "GCTCTTC"
            cdna[sapi_cdna_start0:sapi_cdna_start0 + len(sapi_sequence)] = sapi_sequence
            sapi_genomic_start0 = exon_start0 + sapi_cdna_start0
            chromosome[sapi_genomic_start0:sapi_genomic_start0 + len(sapi_sequence)] = sapi_sequence

        cdna_text = "".join(cdna)
        if strand == -1:
            chromosome[exon_start0:exon_end0] = reverse_complement(cdna_text)
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
            strand=strand,
            cdna=cdna_text,
            cds=cds,
            exons=[
                Exon(
                    start1=exon_start0 + 1,
                    end1=exon_end0,
                    strand=strand,
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

    def test_n_terminal_path_is_sequence_complete(self) -> None:
        client = SyntheticEnsemblClient("human")
        result = design_online(
            species_key="human",
            gene="MockTagGene",
            terminus="N-terminal",
            client=client,
        )
        self.assertTrue(result.sequence_complete)
        self.assertEqual(result.status, "SEQUENCE-COMPLETE COMPUTATIONAL DESIGN")
        self.assertEqual(result.backbone_addgene_id, "169226")
        self.assertEqual(result.donor_payload["payload_length_nt"], 726)
        self.assertEqual(result.donor_payload["tag_length_aa"], 235)
        self.assertEqual(result.donor_payload["linker_peptide"], "GGGGSAS")
        self.assertEqual(result.donor_payload["stop_codon"], "")
        self.assertEqual(result.fusion_protein_length_aa, result.protein_length_aa + 242)
        native_protein = translate(client.record.cds[:-3])
        payload_peptide = str(result.donor_payload["payload_peptide"])
        self.assertEqual(
            result.fusion_protein_sequence,
            native_protein[:1] + payload_peptide + native_protein[1:],
        )
        guide = result.top_guide.target_with_pam
        five = result.five_prime_arm.final_gene_oriented_sequence
        three = result.three_prime_arm.final_gene_oriented_sequence
        self.assertEqual(
            result.cloning_fragments["uha_synthesis_fragment_5to3"],
            NTERM_UHA_PREFIX + guide + five + NTERM_UHA_SUFFIX,
        )
        self.assertEqual(
            result.cloning_fragments["dha_synthesis_fragment_5to3"],
            NTERM_DHA_PREFIX + three + guide + NTERM_DHA_SUFFIX,
        )
        self.assertEqual(
            result.cloning_fragments["expected_sapi_overhangs"],
            {
                "vector_to_uha": "TAC",
                "uha_to_payload": "GTG",
                "payload_to_dha": "AGC",
                "dha_to_vector": "AAT",
            },
        )
        self.assertEqual(result.cloning_fragments["assembled_plasmid_length_nt"], 3947)
        self.assertEqual(result.cloning_fragments["assembled_plasmid_sapi_site_count"], 0)
        coordinates = result.cloning_fragments["assembly_coordinate_map"]
        self.assertEqual((coordinates["donor_start0"], coordinates["donor_end0"]), (1085, 3057))
        self.assertEqual(
            (coordinates["tag_start0"], coordinates["tag_end0"]),
            (1708, 2413),
        )
        self.assertEqual(
            (coordinates["linker_start0"], coordinates["linker_end0"]),
            (2413, 2434),
        )
        plasmid = result.cloning_fragments["assembled_plasmid_5to3"]
        self.assertEqual(plasmid[1082:1085], "TAC")
        self.assertEqual(plasmid[3057:3060], "AAT")
        self.assertEqual(
            result.primer_tail_templates["UHA_reverse_5prime_tail"],
            "CGCGCTCTTCACAC",
        )
        self.assertEqual(
            result.primer_tail_templates["UHA_forward_5prime_tail"],
            "AACGCTCTTCATAC" + guide,
        )
        self.assertEqual(
            result.primer_tail_templates["DHA_forward_5prime_tail"],
            "CGCGCTCTTCGAGC",
        )
        self.assertEqual(
            result.primer_tail_templates["DHA_reverse_5prime_tail"],
            "AACGCTCTTCGATT" + reverse_complement(guide),
        )
        genbank = SeqIO.read(StringIO(assembled_plasmid_genbank(result)), "genbank")
        labels = {
            feature.qualifiers.get("label", [""])[0]
            for feature in genbank.features
        }
        self.assertIn("mNeonGreen", labels)
        self.assertIn("GGGGSAS linker", labels)

    def test_n_terminal_backbone_sequence_and_sapi_architecture(self) -> None:
        self.assertTrue(NTERM_BACKBONE_DNA.is_file())
        metadata = backbone_metadata_for(NTERM_BACKBONE)
        self.assertEqual(metadata["length_nt"], 2765)
        self.assertEqual(
            metadata["snapgene_sha256"],
            "75d9d25b4dac8083c401ee5ac76b080a5f62e42d23357a81b4ac33e84a434177",
        )
        self.assertEqual(
            [site["overhang_5to3"] for site in metadata["sapi_sites"]],
            ["TAC", "GTG", "AGC", "AAT"],
        )
        analysis = n_terminal_backbone_analysis()
        self.assertEqual(
            (
                analysis.vector_left_cut0,
                analysis.payload_start0,
                analysis.payload_end0,
                analysis.vector_right_cut0,
            ),
            (1082, 1117, 1840, 1875),
        )
        self.assertEqual(len(analysis.payload_core_sequence), 723)
        self.assertEqual(len(analysis.payload_sequence), 726)

    def test_reverse_strand_n_terminal_path_is_sequence_complete(self) -> None:
        result = design_online(
            species_key="mouse",
            gene="MockTagGene",
            terminus="N-terminal",
            client=SyntheticEnsemblClient("mouse", strand=-1),
        )
        self.assertTrue(result.sequence_complete)
        self.assertEqual(result.gene_strand, "-")
        self.assertEqual(result.backbone_addgene_id, "169226")
        self.assertEqual(result.cloning_fragments["assembled_plasmid_length_nt"], 3947)

    def test_uploaded_backbone_is_classified_and_assembled_from_sapi_sites(self) -> None:
        custom = infer_custom_backbone_definition(NTERM_BACKBONE_DNA)
        self.assertTrue(custom.is_custom)
        self.assertEqual(custom.terminus, "N-terminal")
        self.assertEqual(custom.payload_length_nt, 726)
        result = design_online(
            species_key="human",
            gene="MockTagGene",
            terminus="N-terminal",
            client=SyntheticEnsemblClient("human"),
            backbone_definition=custom,
        )
        self.assertTrue(result.sequence_complete)
        self.assertEqual(result.backbone_addgene_id, "custom")
        self.assertEqual(
            result.cloning_fragments["expected_sapi_overhangs"],
            custom.overhangs,
        )
        self.assertEqual(result.cloning_fragments["assembled_plasmid_length_nt"], 3947)

        custom_c = infer_custom_backbone_definition(BACKBONE_DNA)
        self.assertEqual(custom_c.terminus, "C-terminal")
        c_result = design_online(
            species_key="human",
            gene="MockTagGene",
            client=SyntheticEnsemblClient("human"),
            backbone_definition=custom_c,
        )
        self.assertTrue(c_result.sequence_complete)
        self.assertEqual(c_result.cloning_fragments["assembled_plasmid_length_nt"], 3950)

    def test_synonymous_pam_blocking_mutation_is_applied(self) -> None:
        result = design_online(
            species_key="human",
            gene="MockTagGene",
            client=SyntheticEnsemblClient("human", guide_requires_blocking=True),
        )
        guide = result.top_guide
        self.assertTrue(result.sequence_complete)
        self.assertFalse(guide.pam_destroyed)
        self.assertTrue(guide.final_pam_destroyed)
        self.assertEqual((guide.pam, guide.final_pam), ("GGG", "GGA"))
        self.assertFalse(guide.blocking_mutation_required)
        blocking = [
            mutation
            for mutation in result.five_prime_arm.mutations
            if mutation.kind == "Guide blocking"
        ]
        self.assertEqual(len(blocking), 1)
        self.assertEqual(
            (
                blocking[0].arm_position1,
                blocking[0].genomic_position1,
                blocking[0].original_codon,
                blocking[0].altered_codon,
                blocking[0].protein_consequence,
            ),
            (600, 2000, "GGG", "GGA", "synonymous (G)"),
        )
        self.assertEqual(
            (blocking[0].pam_before, blocking[0].pam_after),
            ("GGG", "GGA"),
        )
        self.assertEqual(
            result.edited_cds_sequence[result.cds_length_without_stop - 3:result.cds_length_without_stop],
            "GGA",
        )
        coordinates = result.cloning_fragments["assembly_coordinate_map"]
        plasmid = result.cloning_fragments["assembled_plasmid_5to3"]
        self.assertEqual(plasmid[coordinates["uha_start0"] + 599], "A")
        report = design_report(result)
        self.assertIn("PAM before/after: GGG -> GGA", report)
        payload = json.loads(design_json(result))
        self.assertEqual(payload["guides"][0]["final_pam"], "GGA")
        genbank = SeqIO.read(StringIO(assembled_plasmid_genbank(result)), "genbank")
        self.assertIn(
            "Guide blocking",
            {
                feature.qualifiers.get("label", [""])[0]
                for feature in genbank.features
            },
        )

    def test_synonymous_seed_blocking_when_pam_cannot_change(self) -> None:
        result = design_online(
            species_key="human",
            gene="MockTagGene",
            client=SyntheticEnsemblClient(
                "human", guide_requires_seed_blocking=True
            ),
        )
        guide = result.top_guide
        self.assertTrue(result.sequence_complete)
        self.assertEqual((guide.pam, guide.final_pam), ("TGG", "TGG"))
        self.assertFalse(guide.final_pam_destroyed)
        self.assertEqual(guide.final_longest_retained_segment, 13)
        blocking = [
            mutation
            for mutation in result.five_prime_arm.mutations
            if mutation.kind == "Guide blocking"
        ]
        self.assertEqual(len(blocking), 1)
        self.assertEqual(
            (blocking[0].arm_position1, blocking[0].original_codon, blocking[0].altered_codon),
            (591, "CAC", "CAT"),
        )

    def test_no_safe_synonymous_guide_block_stays_blocked(self) -> None:
        result = design_online(
            species_key="mouse",
            gene="MockTagGene",
            client=SyntheticEnsemblClient(
                "mouse", guide_has_no_synonymous_block=True
            ),
        )
        self.assertFalse(result.sequence_complete)
        self.assertEqual(
            result.status,
            "DESIGN BLOCKED - NO SEQUENCE-COMPLETE OUTPUT",
        )
        self.assertTrue(result.top_guide.blocking_mutation_required)
        self.assertIn("No synonymous coding change", result.top_guide.blocking_mutation_note)

    def test_generic_coding_sapi_domestication(self) -> None:
        result = design_online(
            species_key="mouse",
            gene="MockTagGene",
            client=SyntheticEnsemblClient("mouse", coding_sapi_site=True),
        )
        self.assertTrue(result.sequence_complete)
        mutations = [
            mutation
            for mutation in result.five_prime_arm.mutations
            if mutation.kind == "SapI domestication"
        ]
        self.assertEqual(len(mutations), 1)
        self.assertEqual(
            (
                mutations[0].arm_position1,
                mutations[0].genomic_position1,
                mutations[0].original_codon,
                mutations[0].altered_codon,
            ),
            (306, 1706, "CTC", "CTA"),
        )
        self.assertEqual(result.five_prime_arm.final_sapi_sites, [])

    def test_noncoding_sapi_site_is_automatically_domesticated(self) -> None:
        result = design_online(
            species_key="human",
            gene="MockTagGene",
            client=SyntheticEnsemblClient("human", noncoding_sapi_site=True),
        )
        self.assertTrue(result.sequence_complete)
        self.assertEqual(result.three_prime_arm.final_sapi_sites, [])
        mutations = result.three_prime_arm.mutations
        self.assertEqual(len(mutations), 1)
        self.assertEqual(
            (
                mutations[0].kind,
                mutations[0].arm_position1,
                mutations[0].genomic_position1,
                mutations[0].reference_base,
                mutations[0].alternate_base,
                mutations[0].transcript_position1,
                mutations[0].protein_consequence,
            ),
            (
                "SapI domestication",
                104,
                2107,
                "C",
                "A",
                1107,
                "non-coding (not translated)",
            ),
        )
        qc_rows = sapi_qc_rows(result)
        self.assertEqual(len(qc_rows), 1)
        self.assertEqual(qc_rows[0]["Status"], "Resolved")
        self.assertEqual(qc_rows[0]["Codon change"], "Not applicable")
        self.assertEqual(qc_rows[0]["Protein consequence"], "non-coding (not translated)")

    def test_intronic_sapi_fallback_respects_splice_flanks(self) -> None:
        species = SPECIES["mouse"]
        record = TranscriptRecord(
            species=species,
            gene_symbol="MockSpliceGene",
            gene_id="MOCKSPLICE1",
            transcript_id="MOCKSPLICET1",
            transcript_version="1",
            chromosome="1",
            strand=1,
            cdna="A" * 100,
            cds="ATG" + "GCT" * 18 + "TAA",
            exons=[Exon(start1=1001, end1=1100, strand=1, stable_id="MOCKE1")],
            source="synthetic test fixture",
            source_release="test",
        )
        mapping = list(range(1000, 1100))
        site = {"motif": "GCTCTTC", "position0": 0, "position1": 1}
        intronic_arm = HomologyArm(
            name="5-prime homology arm",
            length=7,
            chromosome="1",
            genomic_start0=1200,
            genomic_end0=1207,
            gene_oriented_sequence="GCTCTTC",
            chromosome_forward_sequence="GCTCTTC",
            gc_percent=57.1,
            sapi_sites=[site],
            final_sapi_sites=[site],
        )
        mutated = _automatic_noncoding_sapi_mutation(
            arm=intronic_arm,
            motif="GCTCTTC",
            motif_start0=0,
            record=record,
            transcript_mapping=mapping,
            cds_start_cdna0=10,
            current_cdna=list(record.cdna),
        )
        self.assertIsNotNone(mutated)
        assert mutated is not None
        self.assertEqual(mutated.final_sapi_sites, [])
        self.assertEqual(
            (
                mutated.mutations[0].arm_position1,
                mutated.mutations[0].genomic_position1,
                mutated.mutations[0].reference_base,
                mutated.mutations[0].alternate_base,
                mutated.mutations[0].transcript_position1,
            ),
            (5, 1205, "T", "A", None),
        )

        splice_flank_arm = HomologyArm(
            name="3-prime homology arm",
            length=7,
            chromosome="1",
            genomic_start0=1097,
            genomic_end0=1104,
            gene_oriented_sequence="GCTCTTC",
            chromosome_forward_sequence="GCTCTTC",
            gc_percent=57.1,
            sapi_sites=[site],
            final_sapi_sites=[site],
        )
        self.assertIsNone(
            _automatic_noncoding_sapi_mutation(
                arm=splice_flank_arm,
                motif="GCTCTTC",
                motif_start0=0,
                record=record,
                transcript_mapping=mapping,
                cds_start_cdna0=10,
                current_cdna=list(record.cdna),
            )
        )

    def test_multiple_coding_sapi_sites_are_domesticated(self) -> None:
        result = design_online(
            species_key="human",
            gene="MockTagGene",
            client=SyntheticEnsemblClient(
                "human", multiple_coding_sapi_sites=True
            ),
        )
        self.assertTrue(result.sequence_complete)
        self.assertEqual(result.five_prime_arm.final_sapi_sites, [])
        self.assertEqual(
            len(
                [
                    mutation
                    for mutation in result.five_prime_arm.mutations
                    if mutation.kind == "SapI domestication"
                ]
            ),
            2,
        )
        qc_rows = sapi_qc_rows(result)
        self.assertEqual(len(qc_rows), 2)
        self.assertTrue(all(row["Status"] == "Resolved" for row in qc_rows))

    def test_reverse_strand_guide_blocking_and_sapi_domestication(self) -> None:
        result = design_online(
            species_key="mouse",
            gene="MockTagGene",
            client=SyntheticEnsemblClient(
                "mouse",
                guide_requires_blocking=True,
                coding_sapi_site=True,
                strand=-1,
            ),
        )
        self.assertTrue(result.sequence_complete)
        self.assertEqual(result.gene_strand, "-")
        self.assertEqual(result.top_guide.chromosome_strand, "-")
        self.assertEqual(result.top_guide.final_pam, "GGA")
        self.assertTrue(result.top_guide.final_pam_destroyed)
        self.assertEqual(
            {mutation.kind for mutation in result.five_prime_arm.mutations},
            {"SapI domestication", "Guide blocking"},
        )
        guide_mutation = next(
            mutation
            for mutation in result.five_prime_arm.mutations
            if mutation.kind == "Guide blocking"
        )
        self.assertEqual(
            (guide_mutation.arm_position1, guide_mutation.genomic_position1),
            (600, 1501),
        )


if __name__ == "__main__":
    unittest.main()
