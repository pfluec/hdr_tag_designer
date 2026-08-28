from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from hashlib import sha256
from pathlib import Path
from typing import Any

from .sequence import clean_dna, translate
from .snapgene import SnapGeneDocument, read_snapgene

BACKBONE_NAME = "TVBB C-term-mNeongreen"
BACKBONE_ADDGENE_ID = "169227"
BACKBONE_TERMINUS = "C-terminal"
NTERM_BACKBONE_NAME = "TVBB N-term-mNeongreen"
NTERM_BACKBONE_ADDGENE_ID = "169226"
NTERM_BACKBONE_TERMINUS = "N-terminal"

# Exact C-terminal homology-arm fragment architecture from Bollen et al. S1/S3.
CTERM_UHA_PREFIX = "AACGCTCTTCATAC"
CTERM_UHA_SUFFIX = "GGCTGAAGAGCGCG"
CTERM_DHA_PREFIX = "CGCGCTCTTCGTGA"
CTERM_DHA_SUFFIX = "AATCGAAGAGCGTT"

# PCR-primer 5' tails from the same protocol. Locus-specific annealing sequences
# are intentionally not designed in this prototype.
CTERM_UHA_FORWARD_PRIMER_TAIL_PREFIX = "AACGCTCTTCATAC"
CTERM_UHA_REVERSE_PRIMER_TAIL = "CGCGCTCTTCAGCC"
CTERM_DHA_FORWARD_PRIMER_TAIL = "CGCGCTCTTCGTGA"
CTERM_DHA_REVERSE_PRIMER_TAIL_PREFIX = "AACGCTCTTCGATT"

# Exact N-terminal architecture from Bollen et al. supplementary S1.
NTERM_UHA_PREFIX = "AACGCTCTTCATAC"
NTERM_UHA_SUFFIX = "GTGTGAAGAGCGCG"
NTERM_DHA_PREFIX = "CGCGCTCTTCGAGC"
NTERM_DHA_SUFFIX = "AATCGAAGAGCGTT"
NTERM_UHA_FORWARD_PRIMER_TAIL_PREFIX = "AACGCTCTTCATAC"
NTERM_UHA_REVERSE_PRIMER_TAIL = "CGCGCTCTTCACAC"
NTERM_DHA_FORWARD_PRIMER_TAIL = "CGCGCTCTTCGAGC"
NTERM_DHA_REVERSE_PRIMER_TAIL_PREFIX = "AACGCTCTTCGATT"

SAPI_RECOGNITION_MOTIFS = ("GCTCTTC", "GAAGAGC")
SAPI_OVERHANGS = {
    "vector_to_uha": "TAC",
    "uha_to_payload": "GGC",
    "payload_to_dha": "TGA",
    "dha_to_vector": "AAT",
}
NTERM_SAPI_OVERHANGS = {
    "vector_to_uha": "TAC",
    "uha_to_payload": "GTG",
    "payload_to_dha": "AGC",
    "dha_to_vector": "AAT",
}

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
PAYLOAD_FASTA = DATA_DIR / "bollen_169227_c_terminal_payload.fa"
BACKBONE_DNA = DATA_DIR / "addgene_169227.dna"
NTERM_BACKBONE_DNA = DATA_DIR / "addgene-169226.dna"


@dataclass(frozen=True)
class BackboneDefinition:
    """Sequence architecture for a verified built-in transfer backbone."""

    key: str
    name: str
    addgene_id: str
    terminus: str
    dna_path: Path
    expected_length_nt: int
    expected_sha256: str
    overhang_items: tuple[tuple[str, str], ...]
    uha_prefix: str
    uha_suffix: str
    dha_prefix: str
    dha_suffix: str
    uha_forward_primer_tail_prefix: str
    uha_reverse_primer_tail: str
    dha_forward_primer_tail: str
    dha_reverse_primer_tail_prefix: str
    payload_length_nt: int
    payload_sha256: str
    tag_start0: int
    tag_end0: int
    linker_start0: int
    linker_end0: int
    payload_supplies_stop: bool
    tag_name: str = "mNeonGreen"
    is_custom: bool = False
    fusion_compatible: bool = True
    payload_kind: str = "in-frame fusion"
    payload_warning: str = ""
    source_filename: str = ""

    @property
    def overhangs(self) -> dict[str, str]:
        return dict(self.overhang_items)


CTERM_BACKBONE = BackboneDefinition(
    key="bollen_cterm_mneongreen",
    name=BACKBONE_NAME,
    addgene_id=BACKBONE_ADDGENE_ID,
    terminus=BACKBONE_TERMINUS,
    dna_path=BACKBONE_DNA,
    expected_length_nt=2768,
    expected_sha256="b5ccaa5a257b71a1f2bac05ab15785f07098a46ed805c2862b5beeade04046b1",
    overhang_items=tuple(SAPI_OVERHANGS.items()),
    uha_prefix=CTERM_UHA_PREFIX,
    uha_suffix=CTERM_UHA_SUFFIX,
    dha_prefix=CTERM_DHA_PREFIX,
    dha_suffix=CTERM_DHA_SUFFIX,
    uha_forward_primer_tail_prefix=CTERM_UHA_FORWARD_PRIMER_TAIL_PREFIX,
    uha_reverse_primer_tail=CTERM_UHA_REVERSE_PRIMER_TAIL,
    dha_forward_primer_tail=CTERM_DHA_FORWARD_PRIMER_TAIL,
    dha_reverse_primer_tail_prefix=CTERM_DHA_REVERSE_PRIMER_TAIL_PREFIX,
    payload_length_nt=729,
    payload_sha256="d750c54d2ed969421c84377e7e626edcd902f1eb71ea53fd17b3b12a6e588724",
    tag_start0=21,
    tag_end0=726,
    linker_start0=0,
    linker_end0=21,
    payload_supplies_stop=True,
)

NTERM_BACKBONE = BackboneDefinition(
    key="bollen_nterm_mneongreen",
    name=NTERM_BACKBONE_NAME,
    addgene_id=NTERM_BACKBONE_ADDGENE_ID,
    terminus=NTERM_BACKBONE_TERMINUS,
    dna_path=NTERM_BACKBONE_DNA,
    expected_length_nt=2765,
    expected_sha256="75d9d25b4dac8083c401ee5ac76b080a5f62e42d23357a81b4ac33e84a434177",
    overhang_items=tuple(NTERM_SAPI_OVERHANGS.items()),
    uha_prefix=NTERM_UHA_PREFIX,
    uha_suffix=NTERM_UHA_SUFFIX,
    dha_prefix=NTERM_DHA_PREFIX,
    dha_suffix=NTERM_DHA_SUFFIX,
    uha_forward_primer_tail_prefix=NTERM_UHA_FORWARD_PRIMER_TAIL_PREFIX,
    uha_reverse_primer_tail=NTERM_UHA_REVERSE_PRIMER_TAIL,
    dha_forward_primer_tail=NTERM_DHA_FORWARD_PRIMER_TAIL,
    dha_reverse_primer_tail_prefix=NTERM_DHA_REVERSE_PRIMER_TAIL_PREFIX,
    payload_length_nt=726,
    payload_sha256="8bdcf8abfbcbb281b2148f6e3b509930734d27237ed8564dee651363199646ed",
    tag_start0=0,
    tag_end0=705,
    linker_start0=705,
    linker_end0=726,
    payload_supplies_stop=False,
)


def backbone_for_terminus(terminus: str) -> BackboneDefinition:
    if terminus.upper().startswith("C"):
        return CTERM_BACKBONE
    if terminus.upper().startswith("N"):
        return NTERM_BACKBONE
    raise ValueError("Terminus must be N-terminal or C-terminal")


def infer_custom_backbone_definition(
    path: str | Path,
    *,
    source_filename: str | None = None,
) -> BackboneDefinition:
    """Classify a circular .dna backbone from its four SapI overhangs.

    Custom payloads must retain either the Bollen N- or C-terminal SapI overhang
    architecture. The sequence between the two inner cuts is extracted directly
    from the uploaded file; no payload identity or single-ORF structure is assumed.
    """
    document = read_snapgene(path)
    if document.topology != "circular":
        raise ValueError("Custom transfer backbones must be marked circular")
    sites = find_sapi_sites(document.sequence, circular=True)
    if len(sites) != 4:
        raise ValueError(
            f"A custom transfer backbone must contain exactly four SapI sites; found {len(sites)}"
        )
    observed = [str(site["overhang_5to3"]) for site in sites]
    templates = (CTERM_BACKBONE, NTERM_BACKBONE)
    template = next(
        (
            candidate
            for candidate in templates
            if observed == list(candidate.overhangs.values())
        ),
        None,
    )
    if template is None:
        raise ValueError(
            "Custom backbone SapI overhangs must occur in one of the supported orders: "
            "TAC/GGC/TGA/AAT (C-terminal) or TAC/GTG/AGC/AAT (N-terminal); "
            f"found {'/'.join(observed)}"
        )
    cuts = [int(site["top_strand_cut_boundary0"]) for site in sites]
    if cuts != sorted(cuts):
        raise ValueError(
            "Custom backbone cassette crosses the sequence origin; rotate the circular "
            "sequence so the TAC cut precedes the payload before uploading"
        )
    payload_core = document.sequence[cuts[1]:cuts[2]]
    payload = payload_core + template.overhangs["payload_to_dha"]
    payload_classification = _classify_custom_payload(payload, template)

    uploaded_filename = Path(source_filename or document.path.name).name
    uploaded_stem = Path(uploaded_filename).stem

    return BackboneDefinition(
        key=f"custom_{template.terminus.lower().replace('-', '_')}",
        name=uploaded_stem,
        addgene_id="custom",
        terminus=template.terminus,
        dna_path=document.path,
        expected_length_nt=document.length,
        expected_sha256=document.sha256,
        overhang_items=template.overhang_items,
        uha_prefix=template.uha_prefix,
        uha_suffix=template.uha_suffix,
        dha_prefix=template.dha_prefix,
        dha_suffix=template.dha_suffix,
        uha_forward_primer_tail_prefix=template.uha_forward_primer_tail_prefix,
        uha_reverse_primer_tail=template.uha_reverse_primer_tail,
        dha_forward_primer_tail=template.dha_forward_primer_tail,
        dha_reverse_primer_tail_prefix=template.dha_reverse_primer_tail_prefix,
        payload_length_nt=len(payload),
        payload_sha256=sha256(payload.encode("ascii")).hexdigest(),
        tag_start0=int(payload_classification["tag_start0"]),
        tag_end0=int(payload_classification["tag_end0"]),
        linker_start0=int(payload_classification["linker_start0"]),
        linker_end0=int(payload_classification["linker_end0"]),
        payload_supplies_stop=template.payload_supplies_stop,
        tag_name="custom tag",
        is_custom=True,
        fusion_compatible=bool(payload_classification["fusion_compatible"]),
        payload_kind=str(payload_classification["payload_kind"]),
        payload_warning=str(payload_classification["payload_warning"]),
        source_filename=uploaded_filename,
    )


def _classify_custom_payload(
    payload: str,
    template: BackboneDefinition,
) -> dict[str, int | bool | str]:
    """Describe a custom payload without rejecting complex cassette structures."""
    if template.payload_supplies_stop:
        expected_linker_start0, expected_linker_end0 = 0, 21
        payload_coding = payload[:-3]
    else:
        expected_linker_end0 = len(payload)
        expected_linker_start0 = max(0, expected_linker_end0 - 21)
        payload_coding = payload
    linker_matches = (
        expected_linker_end0 - expected_linker_start0 == 21
        and translate(payload[expected_linker_start0:expected_linker_end0]) == "GGGGSAS"
    )
    frame_divisible = len(payload_coding) % 3 == 0
    translated = translate(payload_coding) if frame_divisible else ""
    unambiguous = "X" not in translated
    no_internal_stop = "*" not in translated
    fusion_compatible = bool(
        linker_matches and frame_divisible and unambiguous and no_internal_stop
    )
    if fusion_compatible and template.payload_supplies_stop:
        linker_start0, linker_end0 = 0, 21
        tag_start0, tag_end0 = 21, len(payload) - 3
    elif fusion_compatible:
        tag_start0, tag_end0 = 0, len(payload) - 21
        linker_start0, linker_end0 = tag_end0, len(payload)
    else:
        linker_start0 = linker_end0 = 0
        tag_start0 = 0
        tag_end0 = len(payload_coding)

    issues: list[str] = []
    if not linker_matches:
        issues.append("the architecture-specific GGGGSAS fusion linker was not detected")
    if not frame_divisible:
        issues.append(
            f"the payload coding span is {len(payload_coding)} nt and is not divisible by three"
        )
    if frame_divisible and not unambiguous:
        issues.append("the payload translation contains ambiguous codons")
    if frame_divisible and not no_internal_stop:
        issues.append("the payload contains one or more in-frame stop codons")
    payload_warning = (
        "Custom payload is retained as a multi-ORF/non-coding or non-frame cassette; "
        + "; ".join(issues)
        + ". No single fusion-protein translation is asserted."
        if issues
        else ""
    )

    return {
        "tag_start0": tag_start0,
        "tag_end0": tag_end0,
        "linker_start0": linker_start0,
        "linker_end0": linker_end0,
        "fusion_compatible": fusion_compatible,
        "payload_kind": "in-frame fusion" if fusion_compatible else "complex cassette",
        "payload_warning": payload_warning,
    }


@dataclass(frozen=True)
class FixedBackboneAnalysis:
    document: SnapGeneDocument
    sapi_sites: tuple[dict[str, Any], ...]
    site_by_overhang: dict[str, dict[str, Any]]
    vector_left_cut0: int
    vector_right_cut0: int
    payload_start0: int
    payload_end0: int
    payload_core_sequence: str
    payload_sequence: str
    payload_matches_supplement: bool


def load_supplemental_c_terminal_payload() -> str:
    lines = PAYLOAD_FASTA.read_text(encoding="utf-8").splitlines()
    payload = clean_dna("".join(line for line in lines if line and not line.startswith(">")))
    if len(payload) != 729:
        raise ValueError(f"Unexpected supplementary #169227 payload length: {len(payload)} nt")
    if payload[:21] != "GGCGGAGGCGGCAGCGCCAGC":
        raise ValueError("Unexpected supplementary #169227 linker sequence")
    if payload[-3:] != "TGA":
        raise ValueError("Unexpected supplementary #169227 stop codon")
    if any(motif in payload for motif in SAPI_RECOGNITION_MOTIFS):
        raise ValueError("The supplementary #169227 payload contains an internal SapI site")
    return payload


def _circular_slice(sequence: str, start0: int, length: int) -> str:
    if not sequence:
        return ""
    size = len(sequence)
    return "".join(sequence[(start0 + offset) % size] for offset in range(length))


def find_sapi_sites(sequence: str, *, circular: bool = False) -> list[dict[str, Any]]:
    """Return SapI recognition sites, top-strand cuts, and 3-nt overhangs.

    SapI is a type-IIS enzyme with the pattern ``GCTCTTC(1/4)``. In a
    chromosome-forward string, a forward ``GCTCTTC`` site cuts the top strand
    eight bases after the first recognition base. A reverse-oriented
    ``GAAGAGC`` site cuts the top strand four bases before the first displayed
    recognition base. Coordinates are 0-based boundaries.
    """
    sequence = clean_dna(sequence)
    size = len(sequence)
    search_sequence = sequence + sequence[:6] if circular and size else sequence
    sites: list[dict[str, Any]] = []
    for motif in SAPI_RECOGNITION_MOTIFS:
        start = 0
        while True:
            position0 = search_sequence.find(motif, start)
            if position0 < 0:
                break
            if position0 >= size:
                break
            if motif == "GCTCTTC":
                orientation = "+"
                top_cut0 = position0 + 8
            else:
                orientation = "-"
                top_cut0 = position0 - 4
            cut_for_slice = top_cut0 % size if circular and size else top_cut0
            if circular:
                overhang = _circular_slice(sequence, cut_for_slice, 3)
            elif 0 <= cut_for_slice <= size - 3:
                overhang = sequence[cut_for_slice:cut_for_slice + 3]
            else:
                overhang = ""
            sites.append(
                {
                    "motif": motif,
                    "orientation": orientation,
                    "recognition_start0": position0,
                    "recognition_end0": position0 + 7,
                    "recognition_interval_1based": f"{position0 + 1}-{position0 + 7}",
                    "top_strand_cut_boundary0": cut_for_slice,
                    "overhang_5to3": overhang,
                }
            )
            start = position0 + 1
    return sorted(sites, key=lambda site: int(site["recognition_start0"]))


@lru_cache(maxsize=4)
def analyze_builtin_backbone(definition: BackboneDefinition) -> FixedBackboneAnalysis:
    document = read_snapgene(definition.dna_path)
    if document.length != definition.expected_length_nt:
        raise ValueError(
            f"Unexpected Addgene #{definition.addgene_id} backbone length: {document.length} nt "
            f"(expected {definition.expected_length_nt})"
        )
    if document.sha256 != definition.expected_sha256:
        raise ValueError(
            f"Addgene #{definition.addgene_id} SnapGene checksum does not match the verified input"
        )
    if document.topology != "circular":
        raise ValueError(
            f"The supplied Addgene #{definition.addgene_id} SnapGene file is not marked circular"
        )

    sites = find_sapi_sites(document.sequence, circular=True)
    if len(sites) != 4:
        raise ValueError(
            f"Expected four SapI sites in Addgene #{definition.addgene_id}, found {len(sites)}"
        )
    site_by_overhang: dict[str, dict[str, Any]] = {}
    for site in sites:
        overhang = str(site["overhang_5to3"])
        if overhang in site_by_overhang:
            raise ValueError(f"SapI overhang {overhang} occurs more than once in the backbone")
        site_by_overhang[overhang] = site
    expected = set(definition.overhangs.values())
    if set(site_by_overhang) != expected:
        raise ValueError(
            f"The uploaded backbone SapI overhangs do not match the Bollen {definition.terminus} "
            f"architecture: found {sorted(site_by_overhang)}, expected {sorted(expected)}"
        )

    overhangs = definition.overhangs
    vector_left_cut0 = int(
        site_by_overhang[overhangs["vector_to_uha"]]["top_strand_cut_boundary0"]
    )
    payload_start0 = int(
        site_by_overhang[overhangs["uha_to_payload"]]["top_strand_cut_boundary0"]
    )
    payload_end0 = int(
        site_by_overhang[overhangs["payload_to_dha"]]["top_strand_cut_boundary0"]
    )
    vector_right_cut0 = int(
        site_by_overhang[overhangs["dha_to_vector"]]["top_strand_cut_boundary0"]
    )
    if not (
        0 < vector_left_cut0 < payload_start0 < payload_end0 < vector_right_cut0 < document.length
    ):
        raise ValueError(
            f"Unexpected SapI cut order in the uploaded Addgene #{definition.addgene_id} backbone"
        )

    # The isolated backbone tag fragment runs from the GGC top-strand cut to the
    # TGA top-strand cut and is 726 nt. In the ligated donor, the TGA cohesive
    # end is supplied once at the payload-to-DHA junction, giving the 729-nt
    # linker-mNeonGreen-stop payload listed in supplementary S2.
    payload_core = document.sequence[payload_start0:payload_end0]
    payload = payload_core + overhangs["payload_to_dha"]
    if len(payload) != definition.payload_length_nt:
        raise ValueError(
            f"Unexpected Addgene #{definition.addgene_id} payload length: {len(payload)} nt "
            f"(expected {definition.payload_length_nt})"
        )
    if sha256(payload.encode("ascii")).hexdigest() != definition.payload_sha256:
        raise ValueError(
            f"The payload reconstructed from Addgene #{definition.addgene_id} does not match "
            "the verified built-in sequence"
        )

    return FixedBackboneAnalysis(
        document=document,
        sapi_sites=tuple(sites),
        site_by_overhang=site_by_overhang,
        vector_left_cut0=vector_left_cut0,
        vector_right_cut0=vector_right_cut0,
        payload_start0=payload_start0,
        payload_end0=payload_end0,
        payload_core_sequence=payload_core,
        payload_sequence=payload,
        payload_matches_supplement=True,
    )


@lru_cache(maxsize=1)
def fixed_backbone_analysis() -> FixedBackboneAnalysis:
    analysis = analyze_builtin_backbone(CTERM_BACKBONE)
    if analysis.payload_sequence != load_supplemental_c_terminal_payload():
        raise ValueError(
            "The linker-mNeonGreen-stop payload reconstructed from the uploaded backbone "
            "does not match the sequence supplied in Bollen supplementary S2."
        )
    return analysis


@lru_cache(maxsize=1)
def n_terminal_backbone_analysis() -> FixedBackboneAnalysis:
    analysis = analyze_builtin_backbone(NTERM_BACKBONE)
    c_terminal_tag = fixed_backbone_analysis().payload_sequence[21:-3]
    payload = analysis.payload_sequence
    if payload[:705] != c_terminal_tag:
        raise ValueError("The N-terminal and C-terminal mNeonGreen tag sequences differ")
    if translate(payload[705:]) != "GGGGSAS":
        raise ValueError("Unexpected Addgene #169226 N-terminal linker sequence")
    if "*" in translate(payload):
        raise ValueError("The Addgene #169226 N-terminal payload contains a stop codon")
    return analysis


def load_c_terminal_payload() -> str:
    """Return the payload extracted directly from the uploaded .dna backbone."""
    return fixed_backbone_analysis().payload_sequence


def load_n_terminal_payload() -> str:
    """Return the tag-linker payload extracted from Addgene #169226."""
    return n_terminal_backbone_analysis().payload_sequence


def _analysis_for_definition(definition: BackboneDefinition) -> FixedBackboneAnalysis:
    if definition == CTERM_BACKBONE:
        return fixed_backbone_analysis()
    if definition == NTERM_BACKBONE:
        return n_terminal_backbone_analysis()
    return analyze_builtin_backbone(definition)


def backbone_metadata_for(definition: BackboneDefinition) -> dict[str, Any]:
    analysis = _analysis_for_definition(definition)
    document = analysis.document
    overhangs = definition.overhangs
    return {
        "key": definition.key,
        "name": definition.name,
        "addgene_id": definition.addgene_id,
        "terminus": definition.terminus,
        "snapgene_file": getattr(definition, "source_filename", "") or document.path.name,
        "snapgene_sha256": document.sha256,
        "length_nt": document.length,
        "topology": document.topology,
        "annotated_feature_count": len(document.features),
        "sapi_site_count": len(analysis.sapi_sites),
        "sapi_sites": [dict(site) for site in analysis.sapi_sites],
        "vector_left_cut_boundary0": analysis.vector_left_cut0,
        "vector_right_cut_boundary0": analysis.vector_right_cut0,
        "payload_core_start0": analysis.payload_start0,
        "payload_core_end0": analysis.payload_end0,
        "payload_core_length_nt": len(analysis.payload_core_sequence),
        "payload_length_nt": len(analysis.payload_sequence),
        "payload_stop_overhang_5to3": overhangs["payload_to_dha"],
        "payload_matches_bollen_s2": (
            analysis.payload_matches_supplement and not definition.is_custom
        ),
        "payload_sequence_verified": analysis.payload_matches_supplement,
        "expected_sapi_overhangs": overhangs,
    }


def backbone_metadata() -> dict[str, Any]:
    """Compatibility wrapper for the original C-terminal backbone API."""
    return backbone_metadata_for(CTERM_BACKBONE)


def payload_metadata_for(
    definition: BackboneDefinition,
) -> dict[str, str | int | bool]:
    analysis = _analysis_for_definition(definition)
    payload = analysis.payload_sequence
    linker_coding = payload[definition.linker_start0:definition.linker_end0]
    tag_coding = payload[definition.tag_start0:definition.tag_end0]
    payload_coding = payload[:-3] if definition.payload_supplies_stop else payload
    architecture = (
        (
            f"GGGGSAS linker + {definition.tag_name} + stop"
            if definition.payload_supplies_stop
            else f"{definition.tag_name} + GGGGSAS linker"
        )
        if definition.fusion_compatible
        else f"Custom {definition.terminus} complex cassette"
    )
    linker_peptide = translate(linker_coding) if linker_coding else ""
    tag_in_frame = definition.fusion_compatible and len(tag_coding) % 3 == 0
    payload_in_frame = definition.fusion_compatible and len(payload_coding) % 3 == 0
    tag_peptide = translate(tag_coding) if tag_in_frame else ""
    payload_peptide = translate(payload_coding) if payload_in_frame else ""
    return {
        "name": architecture,
        "backbone_name": definition.name,
        "addgene_id": definition.addgene_id,
        "terminus": definition.terminus,
        "payload_sequence_5to3": payload,
        "payload_length_nt": len(payload),
        "linker_coding_sequence": linker_coding,
        "linker_peptide": linker_peptide,
        "tag_coding_sequence": tag_coding,
        "tag_name": definition.tag_name,
        "tag_length_nt": len(tag_coding),
        "tag_peptide": tag_peptide,
        "tag_length_aa": len(tag_peptide),
        "payload_coding_sequence": payload_coding,
        "payload_coding_length_nt": len(payload_coding),
        "payload_peptide": payload_peptide,
        "payload_peptide_length_aa": len(payload_peptide),
        "fusion_compatible": definition.fusion_compatible,
        "payload_kind": definition.payload_kind,
        "payload_warning": definition.payload_warning,
        "stop_codon": (
            payload[-3:]
            if definition.payload_supplies_stop and definition.fusion_compatible
            else ""
        ),
        "payload_to_dha_overhang": definition.overhangs["payload_to_dha"],
        "source": (
            f"Extracted from {definition.dna_path.name} and verified against its expected "
            "length, SapI architecture, and SHA-256 digest"
            + (
                "; the conventional fusion frame was also verified"
                if definition.fusion_compatible
                else "; retained as a complex cassette without a single-frame requirement"
            )
        ),
        "matches_bollen_s2": (
            analysis.payload_matches_supplement and not definition.is_custom
        ),
        "sequence_verified": analysis.payload_matches_supplement,
    }


def payload_metadata() -> dict[str, str | int | bool]:
    """Compatibility wrapper for the original C-terminal payload API."""
    return payload_metadata_for(CTERM_BACKBONE)


def _shifted_backbone_features(
    analysis: FixedBackboneAnalysis,
    *,
    definition: BackboneDefinition,
    replacement_length: int,
) -> list[dict[str, Any]]:
    removed_length = analysis.vector_right_cut0 - analysis.vector_left_cut0
    delta = replacement_length - removed_length
    retained: list[dict[str, Any]] = []
    for feature in analysis.document.features:
        start0 = int(feature["start0"])
        end0 = int(feature["end0"])
        if end0 <= analysis.vector_left_cut0:
            shift = 0
        elif start0 >= analysis.vector_right_cut0:
            shift = delta
        else:
            # The original mNeonGreen/Kozak features lie in the replaced cassette.
            continue
        retained.append(
            {
                "type": feature["type"],
                "label": feature["label"],
                "start0": start0 + shift,
                "end0": end0 + shift,
                "strand": int(feature["strand"]),
                "note": (
                    f"Retained from uploaded Addgene #{definition.addgene_id} "
                    "SnapGene annotation"
                ),
            }
        )
    return retained


def assemble_backbone_plasmid(
    definition: BackboneDefinition,
    *,
    target_with_pam: str,
    uha: str,
    dha: str,
) -> dict[str, Any]:
    """Simulate the SapI Golden Gate product in a verified built-in backbone."""
    target = clean_dna(target_with_pam)
    uha = clean_dna(uha)
    dha = clean_dna(dha)
    if len(target) != 23:
        raise ValueError("SpCas9 target-with-PAM must be 23 nt")
    if any(motif in uha or motif in dha for motif in SAPI_RECOGNITION_MOTIFS):
        raise ValueError("Final homology arms must be SapI-free before plasmid assembly")

    analysis = _analysis_for_definition(definition)
    overhangs = definition.overhangs
    payload = analysis.payload_sequence
    donor_insert = target + uha + payload + dha + target
    original = analysis.document.sequence

    # Simulate the top-strand Golden Gate product explicitly. The TAC cohesive
    # end is retained once before the 5-prime donor target, while AAT is already
    # the first three bases of the retained right-hand vector fragment.
    assembled_replacement = overhangs["vector_to_uha"] + donor_insert
    assembled = (
        original[:analysis.vector_left_cut0]
        + assembled_replacement
        + original[analysis.vector_right_cut0:]
    )
    expected_length = (
        len(original)
        - (analysis.vector_right_cut0 - analysis.vector_left_cut0)
        + len(assembled_replacement)
    )
    if len(assembled) != expected_length:
        raise ValueError("Internal plasmid-length calculation failed")

    donor_start0 = analysis.vector_left_cut0 + len(overhangs["vector_to_uha"])
    target_5_start0 = donor_start0
    uha_start0 = target_5_start0 + len(target)
    payload_start0 = uha_start0 + len(uha)
    linker_start0 = payload_start0 + definition.linker_start0
    linker_end0 = payload_start0 + definition.linker_end0
    tag_start0 = payload_start0 + definition.tag_start0
    tag_end0 = payload_start0 + definition.tag_end0
    payload_end0 = payload_start0 + len(payload)
    dha_start0 = payload_end0
    target_3_start0 = dha_start0 + len(dha)
    donor_end0 = target_3_start0 + len(target)

    # The junction overhangs are retained once in the ligated top-strand sequence.
    junctions = {
        "vector_to_uha": {
            "overhang_5to3": overhangs["vector_to_uha"],
            "start0": donor_start0 - 3,
            "end0": donor_start0,
            "observed": assembled[donor_start0 - 3:donor_start0],
            "window_5to3": assembled[donor_start0 - 12:donor_start0 + 12],
        },
        "uha_to_payload": {
            "overhang_5to3": overhangs["uha_to_payload"],
            "start0": payload_start0,
            "end0": payload_start0 + 3,
            "observed": assembled[payload_start0:payload_start0 + 3],
            "window_5to3": assembled[payload_start0 - 12:payload_start0 + 12],
        },
        "payload_to_dha": {
            "overhang_5to3": overhangs["payload_to_dha"],
            "start0": payload_end0 - 3,
            "end0": payload_end0,
            "observed": assembled[payload_end0 - 3:payload_end0],
            "window_5to3": assembled[payload_end0 - 12:payload_end0 + 12],
        },
        "dha_to_vector": {
            "overhang_5to3": overhangs["dha_to_vector"],
            "start0": donor_end0,
            "end0": donor_end0 + 3,
            "observed": assembled[donor_end0:donor_end0 + 3],
            "window_5to3": assembled[donor_end0 - 12:donor_end0 + 12],
        },
    }
    for name, details in junctions.items():
        if details["observed"] != details["overhang_5to3"]:
            raise ValueError(f"Golden Gate junction {name} failed sequence validation")

    final_sapi_sites = find_sapi_sites(assembled, circular=True)
    features = _shifted_backbone_features(
        analysis,
        definition=definition,
        replacement_length=len(assembled_replacement),
    )
    if not definition.fusion_compatible:
        payload_features = [
            {
                "type": "misc_feature",
                "label": "custom payload cassette",
                "start0": payload_start0,
                "end0": payload_end0,
                "strand": 1,
                "note": definition.payload_warning,
            }
        ]
        donor_architecture = "target-UHA-custom-payload-cassette-DHA-target"
        dha_note = "Gene-oriented final arm after custom payload cassette"
    elif definition.payload_supplies_stop:
        payload_features = [
            {
                "type": "misc_feature",
                "label": "GGGGSAS linker",
                "start0": linker_start0,
                "end0": linker_end0,
                "strand": 1,
                "note": "First GGC codon also forms the UHA-to-payload SapI overhang",
            },
            {
                "type": "CDS",
                "label": f"{definition.tag_name} + stop",
                "start0": tag_start0,
                "end0": payload_end0,
                "strand": 1,
                "note": (
                    f"{definition.name} {definition.tag_name} coding sequence "
                    "including TGA stop"
                ),
            },
        ]
        donor_architecture = f"target-UHA-linker-{definition.tag_name}-stop-DHA-target"
        dha_note = "Gene-oriented final arm; endogenous stop omitted"
    else:
        payload_features = [
            {
                "type": "CDS",
                "label": definition.tag_name,
                "start0": tag_start0,
                "end0": tag_end0,
                "strand": 1,
                "note": f"{definition.name} {definition.tag_name} coding sequence",
            },
            {
                "type": "misc_feature",
                "label": "GGGGSAS linker",
                "start0": linker_start0,
                "end0": linker_end0,
                "strand": 1,
                "note": "Final AGC codon also forms the payload-to-DHA SapI overhang",
            },
        ]
        donor_architecture = f"target-UHA-{definition.tag_name}-linker-DHA-target"
        dha_note = "Gene-oriented final arm beginning at endogenous CDS codon 2"
    features.extend(
        [
            {
                "type": "misc_feature",
                "label": "assembled donor cassette",
                "start0": donor_start0,
                "end0": donor_end0,
                "strand": 1,
                "note": donor_architecture,
            },
            {
                "type": "misc_feature",
                "label": "5-prime donor Cas9 target + PAM",
                "start0": target_5_start0,
                "end0": uha_start0,
                "strand": 1,
                "note": target,
            },
            {
                "type": "misc_feature",
                "label": "5-prime homology arm (UHA), final",
                "start0": uha_start0,
                "end0": payload_start0,
                "strand": 1,
                "note": "Gene-oriented final arm including any synonymous corrections",
            },
            {
                "type": "misc_feature",
                "label": "3-prime homology arm (DHA), final",
                "start0": dha_start0,
                "end0": target_3_start0,
                "strand": 1,
                "note": dha_note,
            },
            {
                "type": "misc_feature",
                "label": "3-prime donor Cas9 target + PAM",
                "start0": target_3_start0,
                "end0": donor_end0,
                "strand": 1,
                "note": target,
            },
        ]
    )
    features.extend(payload_features)
    for name, details in junctions.items():
        features.append(
            {
                "type": "misc_feature",
                "label": f"SapI junction {name} ({details['overhang_5to3']})",
                "start0": int(details["start0"]),
                "end0": int(details["end0"]),
                "strand": 1,
                "note": f"Verified 3-nt overhang {details['overhang_5to3']}",
            }
        )
    features.sort(key=lambda feature: (int(feature["start0"]), int(feature["end0"])))

    return {
        "assembled_donor_insert_5to3": donor_insert,
        "assembled_donor_insert_length_nt": len(donor_insert),
        "uploaded_backbone": backbone_metadata_for(definition),
        "backbone_replaced_interval0": {
            "start0": analysis.vector_left_cut0,
            "end0": analysis.vector_right_cut0,
            "length_nt": analysis.vector_right_cut0 - analysis.vector_left_cut0,
        },
        "assembled_replacement_length_nt": len(assembled_replacement),
        "backbone_retained_length_nt": (
            analysis.document.length
            - (analysis.vector_right_cut0 - analysis.vector_left_cut0)
        ),
        "assembled_plasmid_5to3": assembled,
        "assembled_plasmid_length_nt": len(assembled),
        "assembled_plasmid_topology": "circular",
        "assembled_plasmid_sapi_sites": final_sapi_sites,
        "assembled_plasmid_sapi_site_count": len(final_sapi_sites),
        "golden_gate_junctions": junctions,
        "assembled_plasmid_features": features,
        "assembly_coordinate_map": {
            "donor_start0": donor_start0,
            "donor_end0": donor_end0,
            "target_5_start0": target_5_start0,
            "uha_start0": uha_start0,
            "payload_start0": payload_start0,
            "linker_start0": linker_start0,
            "linker_end0": linker_end0,
            "tag_start0": tag_start0,
            "tag_end0": tag_end0,
            "payload_end0": payload_end0,
            "dha_start0": dha_start0,
            "target_3_start0": target_3_start0,
        },
    }


def assemble_c_terminal_plasmid(
    *, target_with_pam: str, uha: str, dha: str
) -> dict[str, Any]:
    """Compatibility wrapper for Addgene #169227 assembly."""
    return assemble_backbone_plasmid(
        CTERM_BACKBONE,
        target_with_pam=target_with_pam,
        uha=uha,
        dha=dha,
    )


def assemble_n_terminal_plasmid(
    *, target_with_pam: str, uha: str, dha: str
) -> dict[str, Any]:
    return assemble_backbone_plasmid(
        NTERM_BACKBONE,
        target_with_pam=target_with_pam,
        uha=uha,
        dha=dha,
    )


def synthesis_fragments_for_backbone(
    definition: BackboneDefinition,
    *,
    target_with_pam: str,
    uha: str,
    dha: str,
) -> dict[str, Any]:
    target = clean_dna(target_with_pam)
    uha = clean_dna(uha)
    dha = clean_dna(dha)
    if len(target) != 23:
        raise ValueError("SpCas9 target-with-PAM must be 23 nt")
    uha_fragment = definition.uha_prefix + target + uha + definition.uha_suffix
    dha_fragment = definition.dha_prefix + dha + target + definition.dha_suffix
    uha_sapi_sites = find_sapi_sites(uha_fragment)
    dha_sapi_sites = find_sapi_sites(dha_fragment)
    overhangs = definition.overhangs
    expected_uha = [overhangs["vector_to_uha"], overhangs["uha_to_payload"]]
    expected_dha = [overhangs["payload_to_dha"], overhangs["dha_to_vector"]]
    if [site["overhang_5to3"] for site in uha_sapi_sites] != expected_uha:
        raise ValueError(
            "The UHA synthesis fragment does not generate the expected "
            f"{'/'.join(expected_uha)} SapI overhangs"
        )
    if [site["overhang_5to3"] for site in dha_sapi_sites] != expected_dha:
        raise ValueError(
            "The DHA synthesis fragment does not generate the expected "
            f"{'/'.join(expected_dha)} SapI overhangs"
        )
    assembly = assemble_backbone_plasmid(
        definition,
        target_with_pam=target,
        uha=uha,
        dha=dha,
    )
    return {
        "uha_synthesis_fragment_5to3": uha_fragment,
        "uha_synthesis_fragment_length_nt": len(uha_fragment),
        "dha_synthesis_fragment_5to3": dha_fragment,
        "dha_synthesis_fragment_length_nt": len(dha_fragment),
        "synthesis_fragment_sapi_sites": {
            "UHA": uha_sapi_sites,
            "DHA": dha_sapi_sites,
        },
        "expected_sapi_overhangs": overhangs,
        **assembly,
    }


def c_terminal_synthesis_fragments(
    *, target_with_pam: str, uha: str, dha: str
) -> dict[str, Any]:
    return synthesis_fragments_for_backbone(
        CTERM_BACKBONE,
        target_with_pam=target_with_pam,
        uha=uha,
        dha=dha,
    )


def n_terminal_synthesis_fragments(
    *, target_with_pam: str, uha: str, dha: str
) -> dict[str, Any]:
    return synthesis_fragments_for_backbone(
        NTERM_BACKBONE,
        target_with_pam=target_with_pam,
        uha=uha,
        dha=dha,
    )
