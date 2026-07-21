from __future__ import annotations

from dataclasses import asdict, dataclass, field
from typing import Any


@dataclass(frozen=True)
class SpeciesConfig:
    key: str
    label: str
    ensembl_name: str
    assembly: str


@dataclass(frozen=True)
class Exon:
    start1: int
    end1: int
    strand: int
    stable_id: str = ""

    @property
    def length(self) -> int:
        return self.end1 - self.start1 + 1


@dataclass
class TranscriptRecord:
    species: SpeciesConfig
    gene_symbol: str
    gene_id: str
    transcript_id: str
    transcript_version: str
    chromosome: str
    strand: int
    cdna: str
    cds: str
    exons: list[Exon]
    source: str
    source_release: str = ""
    notes: list[str] = field(default_factory=list)

    @property
    def display_transcript_id(self) -> str:
        return (
            f"{self.transcript_id}.{self.transcript_version}"
            if self.transcript_version
            else self.transcript_id
        )


@dataclass
class GuideCandidate:
    rank: int
    spacer: str
    pam: str
    chromosome_strand: str
    target_start0: int
    target_end0: int
    pam_start0: int
    pam_end0: int
    nick_boundary0: int
    distance_to_insertion: int
    gc_percent: float
    poly_t: bool
    insertion_disrupts_target: bool
    deletion_disrupts_target: bool
    target_destroyed: bool
    pam_destroyed: bool
    longest_retained_segment: int
    final_longest_retained_segment: int
    blocking_mutation_required: bool
    activity_heuristic: str
    rationale: str
    blocking_mutation_note: str = ""
    final_pam: str = ""
    final_pam_destroyed: bool = False

    @property
    def target_with_pam(self) -> str:
        return f"{self.spacer}{self.pam}"

    @property
    def target_interval_1based(self) -> str:
        return f"{self.target_start0 + 1}-{self.target_end0}"


@dataclass
class SequenceMutation:
    kind: str
    arm_name: str
    arm_position1: int
    genomic_position1: int
    transcript_position1: int | None
    reference_base: str
    alternate_base: str
    original_codon: str = ""
    altered_codon: str = ""
    amino_acid: str = ""
    protein_consequence: str = ""
    pam_before: str = ""
    pam_after: str = ""
    longest_retained_before: int | None = None
    longest_retained_after: int | None = None
    automatic: bool = True
    reason: str = ""


@dataclass
class HomologyArm:
    name: str
    length: int
    chromosome: str
    genomic_start0: int
    genomic_end0: int
    gene_oriented_sequence: str
    chromosome_forward_sequence: str
    gc_percent: float
    sapi_sites: list[dict[str, Any]] = field(default_factory=list)
    corrected_gene_oriented_sequence: str = ""
    final_sapi_sites: list[dict[str, Any]] = field(default_factory=list)
    mutations: list[SequenceMutation] = field(default_factory=list)
    correction_note: str = ""

    @property
    def genomic_interval_1based(self) -> str:
        return f"{self.chromosome}:{self.genomic_start0 + 1}-{self.genomic_end0}"

    @property
    def final_gene_oriented_sequence(self) -> str:
        return self.corrected_gene_oriented_sequence or self.gene_oriented_sequence


@dataclass
class DesignResult:
    status: str
    sequence_complete: bool
    species_label: str
    assembly: str
    gene_symbol: str
    gene_id: str
    transcript_id: str
    chromosome: str
    gene_strand: str
    terminus: str
    nuclease_mode: str
    backbone_name: str
    backbone_addgene_id: str
    insertion_boundary0: int
    removed_genomic_interval_start0: int
    removed_genomic_interval_end0: int
    removed_sequence_gene_oriented: str
    cds_length_without_stop: int
    protein_length_aa: int
    guide_search_radius: int
    homology_arm_length: int
    guide_safety_cutoff_nt: int
    guide_scoring_note: str
    guides: list[GuideCandidate]
    five_prime_arm: HomologyArm
    three_prime_arm: HomologyArm
    donor_payload: dict[str, Any]
    cloning_fragments: dict[str, Any]
    primer_tail_templates: dict[str, str]
    edited_cds_sequence: str
    fusion_protein_sequence: str
    fusion_protein_length_aa: int
    junctions: dict[str, str]
    validations: list[dict[str, str]]
    warnings: list[str]
    provenance: list[str]
    custom_backbones_supported: bool = False

    @property
    def top_guide(self) -> GuideCandidate:
        if not self.guides:
            raise ValueError("No guides were found")
        return self.guides[0]

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)
