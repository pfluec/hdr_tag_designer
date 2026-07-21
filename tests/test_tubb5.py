from __future__ import annotations

import json
from io import StringIO
import unittest

from Bio import SeqIO

from hdr_designer.backbones import (
    BACKBONE_DNA,
    CTERM_DHA_PREFIX,
    CTERM_DHA_SUFFIX,
    CTERM_UHA_PREFIX,
    CTERM_UHA_SUFFIX,
    SAPI_RECOGNITION_MOTIFS,
    backbone_metadata,
    find_sapi_sites,
    fixed_backbone_analysis,
)
from hdr_designer.design import design_tubb5_fixture
from hdr_designer.exports import (
    arms_fasta,
    assembled_plasmid_genbank,
    design_json,
    design_report,
    genotyping_primers_csv,
    locus_context_genbank,
    guides_csv,
    sapi_qc_rows,
)
from hdr_designer.sequence import reverse_complement, translate


class Tubb5DesignTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.result = design_tubb5_fixture()

    def test_sequence_complete_status(self) -> None:
        self.assertTrue(self.result.sequence_complete)
        self.assertEqual(self.result.status, "SEQUENCE-COMPLETE COMPUTATIONAL DESIGN")
        self.assertTrue(
            all(item["status"] in {"PASS", "WARNING", "N/A"} for item in self.result.validations)
        )

    def test_selected_guide(self) -> None:
        guide = self.result.top_guide
        self.assertEqual(guide.spacer, "GAGGCAGAAGAGGAGGCCTA")
        self.assertEqual(guide.pam, "AGG")
        self.assertEqual(guide.target_with_pam, "GAGGCAGAAGAGGAGGCCTAAGG")
        self.assertEqual(guide.chromosome_strand, "-")
        self.assertEqual((guide.target_start0, guide.target_end0), (36_145_871, 36_145_894))
        self.assertEqual(guide.nick_boundary0, 36_145_877)
        self.assertEqual(guide.distance_to_insertion, 1)
        self.assertTrue(guide.target_destroyed)
        self.assertTrue(guide.pam_destroyed)
        self.assertEqual(guide.longest_retained_segment, 18)
        self.assertEqual(guide.final_longest_retained_segment, 18)
        self.assertFalse(guide.blocking_mutation_required)
        self.assertIn("No extra guide-blocking mutation is required", guide.blocking_mutation_note)
        self.assertEqual(
            guide.final_target_with_pam_after_point_mutations,
            guide.target_with_pam,
        )
        self.assertEqual(
            guide.edited_target_region_display,
            "GAGGCAGAAGAGGAGGCC[INSERT 729 nt]GG",
        )
        self.assertEqual(guide.edited_target_deleted_bases, "TAA")
        self.assertEqual(len(guide.edited_target_region_5to3), 23 - 3 + 729)

    def test_genotyping_primer_geometry(self) -> None:
        primer_result = self.result.genotyping_primers
        self.assertEqual(primer_result["status"], "PASS")
        self.assertEqual(primer_result["assembly"], "GRCm39")
        self.assertEqual(set(primer_result["assays"]), {
            "wild_type_locus",
            "five_prime_junction",
            "three_prime_junction",
        })
        self.assertTrue(
            all(assay["status"] == "PASS" for assay in primer_result["assays"].values())
        )

        five_assay = primer_result["assays"]["five_prime_junction"]
        five_genomic = five_assay["forward_primer"]
        five_payload = five_assay["reverse_primer"]
        self.assertTrue(five_genomic["outside_homology_arm"])
        self.assertFalse(five_genomic["present_in_donor_plasmid"])
        self.assertGreaterEqual(
            five_payload["distance_from_5prime_junction_nt"], 150
        )
        self.assertTrue(five_payload["reusable_payload_primer"])

        three_assay = primer_result["assays"]["three_prime_junction"]
        three_payload = three_assay["forward_primer"]
        three_genomic = three_assay["reverse_primer"]
        self.assertGreaterEqual(
            three_payload["distance_from_3prime_junction_nt"], 150
        )
        self.assertTrue(three_genomic["outside_homology_arm"])
        self.assertFalse(three_genomic["present_in_donor_plasmid"])

        wt_assay = primer_result["assays"]["wild_type_locus"]
        self.assertTrue(wt_assay["forward_primer"]["outside_homology_arm"])
        self.assertTrue(wt_assay["reverse_primer"]["outside_homology_arm"])
        self.assertEqual(
            wt_assay["expected_edited_product_size_bp"]
            - wt_assay["expected_wild_type_product_size_bp"],
            729 - 3,
        )
        self.assertEqual(
            len(wt_assay["expected_wild_type_amplicon_sequence_5to3"]),
            wt_assay["expected_wild_type_product_size_bp"],
        )
        self.assertEqual(
            len(wt_assay["expected_edited_amplicon_sequence_5to3"]),
            wt_assay["expected_edited_product_size_bp"],
        )
        self.assertIn(
            self.result.five_prime_arm.gene_oriented_sequence,
            wt_assay["expected_wild_type_amplicon_sequence_5to3"],
        )
        self.assertNotIn(
            self.result.five_prime_arm.final_gene_oriented_sequence,
            wt_assay["expected_wild_type_amplicon_sequence_5to3"],
        )
        self.assertIn(
            self.result.five_prime_arm.final_gene_oriented_sequence,
            wt_assay["expected_edited_amplicon_sequence_5to3"],
        )
        plasmid = self.result.cloning_fragments["assembled_plasmid_5to3"]
        for assay in primer_result["assays"].values():
            for primer in (assay["forward_primer"], assay["reverse_primer"]):
                if primer.get("outside_homology_arm"):
                    self.assertNotIn(primer["sequence_5to3"], plasmid)
                    self.assertNotIn(reverse_complement(primer["sequence_5to3"]), plasmid)

    def test_annotated_wild_type_and_edited_locus_contexts(self) -> None:
        contexts = self.result.locus_contexts
        self.assertEqual(set(contexts), {"orientation", "external_flank_length_nt", "wild_type", "edited"})
        self.assertEqual(contexts["external_flank_length_nt"], 300)
        wild_type = contexts["wild_type"]
        edited = contexts["edited"]
        self.assertEqual(wild_type["length_nt"], 300 + 600 + 3 + 600 + 300)
        self.assertEqual(edited["length_nt"], 300 + 600 + 729 + 600 + 300)
        self.assertEqual(len(wild_type["sequence_5to3"]), wild_type["length_nt"])
        self.assertEqual(len(edited["sequence_5to3"]), edited["length_nt"])
        self.assertIn(
            self.result.five_prime_arm.gene_oriented_sequence,
            wild_type["sequence_5to3"],
        )
        self.assertIn(
            self.result.five_prime_arm.final_gene_oriented_sequence,
            edited["sequence_5to3"],
        )
        self.assertIn(
            self.result.donor_payload["payload_sequence_5to3"],
            edited["sequence_5to3"],
        )
        wt_labels = {item["label"] for item in wild_type["features"]}
        edited_labels = {item["label"] for item in edited["features"]}
        self.assertIn("5-prime homology arm (UHA), reference", wt_labels)
        self.assertIn("3-prime homology arm (DHA), reference", wt_labels)
        self.assertIn("inserted payload", edited_labels)
        self.assertIn("five_prime_junction reverse primer", edited_labels)
        self.assertNotIn("five_prime_junction reverse primer", wt_labels)

        for context_name, expected_length in (
            ("wild_type", wild_type["length_nt"]),
            ("edited", edited["length_nt"]),
        ):
            record = SeqIO.read(StringIO(locus_context_genbank(self.result, context_name)), "genbank")
            self.assertEqual(len(record.seq), expected_length)
            self.assertEqual(record.annotations["topology"], "linear")
            primer_features = [feature for feature in record.features if feature.type == "primer_bind"]
            self.assertTrue(primer_features)
            for feature in primer_features:
                note = feature.qualifiers.get("note", [""])[0]
                primer_sequence = note.split("5'-", 1)[1].split("-3'", 1)[0]
                self.assertEqual(str(feature.extract(record.seq)), primer_sequence)

    def test_insertion_stop_and_arms(self) -> None:
        self.assertEqual(self.result.insertion_boundary0, 36_145_876)
        self.assertEqual(
            (self.result.removed_genomic_interval_start0, self.result.removed_genomic_interval_end0),
            (36_145_873, 36_145_876),
        )
        self.assertEqual(self.result.removed_sequence_gene_oriented, "TAA")
        five = self.result.five_prime_arm
        three = self.result.three_prime_arm
        self.assertEqual(five.length, 600)
        self.assertEqual(three.length, 600)
        self.assertEqual((five.genomic_start0, five.genomic_end0), (36_145_876, 36_146_476))
        self.assertEqual((three.genomic_start0, three.genomic_end0), (36_145_273, 36_145_873))
        self.assertEqual(reverse_complement(five.gene_oriented_sequence), five.chromosome_forward_sequence)
        self.assertEqual(reverse_complement(three.gene_oriented_sequence), three.chromosome_forward_sequence)

    def test_one_silent_sapi_domestication_mutation(self) -> None:
        five = self.result.five_prime_arm
        self.assertEqual(five.sapi_sites, [{"motif": "GCTCTTC", "position0": 395, "position1": 396}])
        self.assertEqual(five.final_sapi_sites, [])
        self.assertEqual(len(five.mutations), 1)
        self.assertEqual(
            [(m.kind, m.arm_position1, m.genomic_position1, m.original_codon, m.altered_codon)
             for m in five.mutations],
            [
                ("SapI domestication", 396, 36_146_081, "GAG", "GAA"),
            ],
        )
        for mutation in five.mutations:
            self.assertEqual(translate(mutation.original_codon), translate(mutation.altered_codon))
        for motif in SAPI_RECOGNITION_MOTIFS:
            self.assertNotIn(motif, five.final_gene_oriented_sequence)
            self.assertNotIn(motif, self.result.three_prime_arm.final_gene_oriented_sequence)

    def test_sapi_quality_control_summary(self) -> None:
        rows = sapi_qc_rows(self.result)
        self.assertEqual(len(rows), 1)
        self.assertEqual(
            {
                "Arm": rows[0]["Arm"],
                "Motif": rows[0]["Motif"],
                "Arm interval": rows[0]["Arm interval"],
                "Genomic interval": rows[0]["Genomic interval"],
                "Status": rows[0]["Status"],
                "Codon change": rows[0]["Codon change"],
                "Protein consequence": rows[0]["Protein consequence"],
            },
            {
                "Arm": "5-prime homology arm",
                "Motif": "GCTCTTC",
                "Arm interval": "396-402",
                "Genomic interval": "chr17:36,146,075-36,146,081",
                "Status": "Resolved",
                "Codon change": "GAG>GAA",
                "Protein consequence": "Glu (E)",
            },
        )
        self.assertIn("arm base 396", rows[0]["Mutation(s)"])
        report = design_report(self.result)
        self.assertIn("SAPI ARM QUALITY CONTROL", report)
        self.assertIn("Original SapI sites found in both arms: 1", report)
        self.assertIn("Original SapI sites resolved: 1", report)

    def test_fixed_payload_and_fusion(self) -> None:
        payload = self.result.donor_payload
        self.assertEqual(payload["payload_length_nt"], 729)
        self.assertEqual(payload["linker_peptide"], "GGGGSAS")
        self.assertEqual(payload["tag_length_nt"], 705)
        self.assertEqual(payload["tag_length_aa"], 235)
        self.assertEqual(payload["stop_codon"], "TGA")
        self.assertEqual(self.result.cds_length_without_stop, 1332)
        self.assertEqual(self.result.protein_length_aa, 444)
        self.assertEqual(len(self.result.edited_cds_sequence), 2058)
        self.assertEqual(self.result.fusion_protein_length_aa, 686)
        self.assertNotIn("*", self.result.fusion_protein_sequence)

    def test_exact_bollen_fragment_architecture(self) -> None:
        guide = self.result.top_guide.target_with_pam
        five = self.result.five_prime_arm.final_gene_oriented_sequence
        three = self.result.three_prime_arm.final_gene_oriented_sequence
        uha = self.result.cloning_fragments["uha_synthesis_fragment_5to3"]
        dha = self.result.cloning_fragments["dha_synthesis_fragment_5to3"]
        donor = self.result.cloning_fragments["assembled_donor_insert_5to3"]
        payload = self.result.donor_payload["payload_sequence_5to3"]
        self.assertEqual(uha, CTERM_UHA_PREFIX + guide + five + CTERM_UHA_SUFFIX)
        self.assertEqual(dha, CTERM_DHA_PREFIX + three + guide + CTERM_DHA_SUFFIX)
        self.assertEqual(len(uha), 651)
        self.assertEqual(len(dha), 651)
        self.assertEqual(donor, guide + five + payload + three + guide)
        self.assertEqual(len(donor), 1975)
        self.assertEqual(
            self.result.cloning_fragments["expected_sapi_overhangs"],
            {
                "vector_to_uha": "TAC",
                "uha_to_payload": "GGC",
                "payload_to_dha": "TGA",
                "dha_to_vector": "AAT",
            },
        )
        self.assertEqual(
            {
                name: [site["overhang_5to3"] for site in sites]
                for name, sites in self.result.cloning_fragments[
                    "synthesis_fragment_sapi_sites"
                ].items()
            },
            {"UHA": ["TAC", "GGC"], "DHA": ["TGA", "AAT"]},
        )
        self.assertEqual(
            [site["overhang_5to3"] for site in find_sapi_sites(uha)],
            ["TAC", "GGC"],
        )
        self.assertEqual(
            [site["overhang_5to3"] for site in find_sapi_sites(dha)],
            ["TGA", "AAT"],
        )

    def test_uploaded_snapgene_backbone(self) -> None:
        self.assertTrue(BACKBONE_DNA.is_file())
        metadata = backbone_metadata()
        self.assertEqual(metadata["length_nt"], 2768)
        self.assertEqual(metadata["topology"], "circular")
        self.assertEqual(metadata["sapi_site_count"], 4)
        self.assertEqual(
            metadata["snapgene_sha256"],
            "b5ccaa5a257b71a1f2bac05ab15785f07098a46ed805c2862b5beeade04046b1",
        )
        self.assertTrue(metadata["payload_matches_bollen_s2"])
        self.assertEqual(
            [site["overhang_5to3"] for site in metadata["sapi_sites"]],
            ["TAC", "GGC", "TGA", "AAT"],
        )
        analysis = fixed_backbone_analysis()
        self.assertEqual(
            (
                analysis.vector_left_cut0,
                analysis.payload_start0,
                analysis.payload_end0,
                analysis.vector_right_cut0,
            ),
            (1082, 1117, 1843, 1878),
        )
        self.assertEqual(len(analysis.payload_core_sequence), 726)
        self.assertEqual(len(analysis.payload_sequence), 729)

    def test_full_circular_plasmid_assembly(self) -> None:
        fragments = self.result.cloning_fragments
        plasmid = fragments["assembled_plasmid_5to3"]
        self.assertEqual(len(plasmid), 3950)
        self.assertEqual(fragments["assembled_plasmid_length_nt"], 3950)
        self.assertEqual(fragments["assembled_plasmid_topology"], "circular")
        self.assertEqual(fragments["assembled_plasmid_sapi_site_count"], 0)
        self.assertEqual(fragments["assembled_plasmid_sapi_sites"], [])
        self.assertEqual(
            {name: value["observed"] for name, value in fragments["golden_gate_junctions"].items()},
            {
                "vector_to_uha": "TAC",
                "uha_to_payload": "GGC",
                "payload_to_dha": "TGA",
                "dha_to_vector": "AAT",
            },
        )
        coordinates = fragments["assembly_coordinate_map"]
        self.assertEqual(coordinates["donor_start0"], 1085)
        self.assertEqual(coordinates["donor_end0"], 3060)
        self.assertEqual(plasmid[1082:1085], "TAC")
        self.assertEqual(plasmid[3060:3063], "AAT")

    def test_genbank_round_trip(self) -> None:
        text = assembled_plasmid_genbank(self.result)
        record = SeqIO.read(StringIO(text), "genbank")
        self.assertEqual(len(record.seq), 3950)
        self.assertEqual(record.annotations["topology"], "circular")
        labels = {
            feature.qualifiers.get("label", [""])[0]
            for feature in record.features
        }
        self.assertIn("5-prime homology arm (UHA), final", labels)
        self.assertIn("mNeonGreen + stop", labels)
        self.assertIn("3-prime homology arm (DHA), final", labels)
        self.assertIn("SapI domestication", labels)
        self.assertNotIn("Guide-blocking mutation", labels)

    def test_selected_target_is_absent_from_edited_locus(self) -> None:
        edited = self.result.junctions["edited_locus_window_5to3"]
        target = self.result.top_guide.target_with_pam
        self.assertNotIn(target, edited)
        self.assertNotIn(reverse_complement(target), edited)

    def test_primer_tail_templates(self) -> None:
        tails = self.result.primer_tail_templates
        self.assertEqual(
            tails["UHA_forward_5prime_tail"],
            "AACGCTCTTCATACGAGGCAGAAGAGGAGGCCTAAGG",
        )
        self.assertEqual(tails["UHA_reverse_5prime_tail"], "CGCGCTCTTCAGCC")
        self.assertEqual(tails["DHA_forward_5prime_tail"], "CGCGCTCTTCGTGA")
        self.assertEqual(
            tails["DHA_reverse_5prime_tail"],
            "AACGCTCTTCGATTCCTTAGGCCTCCTCTTCTGCCTC",
        )
        cloning = self.result.cloning_primers
        self.assertEqual(set(cloning["primers"]), {
            "UHA_forward", "UHA_reverse", "DHA_forward", "DHA_reverse"
        })
        for primer in cloning["primers"].values():
            self.assertEqual(
                primer["full_sequence_5to3"],
                primer["tail_sequence_5to3"] + primer["annealing_sequence_5to3"],
            )
            self.assertGreaterEqual(primer["annealing_length_nt"], 18)
        self.assertTrue(
            cloning["primers"]["UHA_forward"]["annealing_sequence_5to3"]
            == self.result.five_prime_arm.final_gene_oriented_sequence[:20]
        )
        self.assertEqual(
            cloning["primers"]["UHA_reverse"]["annealing_sequence_5to3"],
            reverse_complement(self.result.five_prime_arm.final_gene_oriented_sequence[-20:]),
        )
        self.assertIn("internal mutation", " ".join(cloning["warnings"]))

    def test_exports(self) -> None:
        report = design_report(self.result)
        self.assertIn("SEQUENCE-COMPLETE COMPUTATIONAL DESIGN", report)
        self.assertIn("No extra guide-blocking mutation is required", report)
        self.assertNotIn("Guide-blocking mutation: arm base", report)
        self.assertIn("uha_synthesis_fragment_5to3", report)
        self.assertIn("Simulated final plasmid: 3950 bp", report)
        self.assertIn("GAGGCAGAAGAGGAGGCCTA", guides_csv(self.result))
        primer_csv = genotyping_primers_csv(self.result)
        self.assertIn("five_prime_junction_forward", primer_csv)
        self.assertIn("three_prime_junction_reverse", primer_csv)
        self.assertIn("GENOTYPING PCR ASSAYS", report)
        fasta = arms_fasta(self.result)
        self.assertIn(">Tubb5_5_prime_homology_arm_gene_oriented_FINAL", fasta)
        self.assertIn(">Tubb5_UHA_synthesis_fragment_FINAL", fasta)
        self.assertIn("assembled_circular_plasmid", fasta)
        self.assertIn("five_prime_junction_expected_amplicon", fasta)
        self.assertIn("wild_type_locus_with_300bp_external_flanks", fasta)
        self.assertIn("edited_locus_with_300bp_external_flanks", fasta)
        payload = json.loads(design_json(self.result))
        self.assertTrue(payload["sequence_complete"])
        self.assertEqual(payload["guides"][0]["final_longest_retained_segment"], 18)
        self.assertFalse(payload["guides"][0]["blocking_mutation_required"])
        self.assertEqual(payload["genotyping_primers"]["status"], "PASS")


if __name__ == "__main__":
    unittest.main()
