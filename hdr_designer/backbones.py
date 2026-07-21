from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
from pathlib import Path
from typing import Any

from .sequence import clean_dna, translate
from .snapgene import SnapGeneDocument, read_snapgene

BACKBONE_NAME = "TVBB C-term-mNeongreen"
BACKBONE_ADDGENE_ID = "169227"
BACKBONE_TERMINUS = "C-terminal"

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

# N-terminal arm adapters are retained for locus-preview mode. The fixed
# Addgene #169227 payload is C-terminal and is not released for N-terminal use.
NTERM_UHA_PREFIX = "AACGCTCTTCATAC"
NTERM_UHA_SUFFIX = "GTGTGAAGAGCGCG"
NTERM_DHA_PREFIX = "CGCGCTCTTCGAGC"
NTERM_DHA_SUFFIX = "AATCGAAGAGCGTT"

SAPI_RECOGNITION_MOTIFS = ("GCTCTTC", "GAAGAGC")
SAPI_OVERHANGS = {
    "vector_to_uha": "TAC",
    "uha_to_payload": "GGC",
    "payload_to_dha": "TGA",
    "dha_to_vector": "AAT",
}

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
PAYLOAD_FASTA = DATA_DIR / "bollen_169227_c_terminal_payload.fa"
BACKBONE_DNA = DATA_DIR / "addgene_169227.dna"


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


@lru_cache(maxsize=1)
def fixed_backbone_analysis() -> FixedBackboneAnalysis:
    document = read_snapgene(BACKBONE_DNA)
    if document.length != 2768:
        raise ValueError(
            f"Unexpected Addgene #169227 backbone length: {document.length} nt (expected 2768)"
        )
    if document.topology != "circular":
        raise ValueError("The supplied Addgene #169227 SnapGene file is not marked circular")

    sites = find_sapi_sites(document.sequence, circular=True)
    if len(sites) != 4:
        raise ValueError(
            f"Expected four SapI sites in Addgene #169227, found {len(sites)}"
        )
    site_by_overhang: dict[str, dict[str, Any]] = {}
    for site in sites:
        overhang = str(site["overhang_5to3"])
        if overhang in site_by_overhang:
            raise ValueError(f"SapI overhang {overhang} occurs more than once in the backbone")
        site_by_overhang[overhang] = site
    expected = set(SAPI_OVERHANGS.values())
    if set(site_by_overhang) != expected:
        raise ValueError(
            "The uploaded backbone SapI overhangs do not match the Bollen C-terminal "
            f"architecture: found {sorted(site_by_overhang)}, expected {sorted(expected)}"
        )

    vector_left_cut0 = int(site_by_overhang["TAC"]["top_strand_cut_boundary0"])
    payload_start0 = int(site_by_overhang["GGC"]["top_strand_cut_boundary0"])
    payload_end0 = int(site_by_overhang["TGA"]["top_strand_cut_boundary0"])
    vector_right_cut0 = int(site_by_overhang["AAT"]["top_strand_cut_boundary0"])
    if not (
        0 < vector_left_cut0 < payload_start0 < payload_end0 < vector_right_cut0 < document.length
    ):
        raise ValueError("Unexpected SapI cut order in the uploaded Addgene #169227 backbone")

    # The isolated backbone tag fragment runs from the GGC top-strand cut to the
    # TGA top-strand cut and is 726 nt. In the ligated donor, the TGA cohesive
    # end is supplied once at the payload-to-DHA junction, giving the 729-nt
    # linker-mNeonGreen-stop payload listed in supplementary S2.
    payload_core = document.sequence[payload_start0:payload_end0]
    supplemental = load_supplemental_c_terminal_payload()
    payload = payload_core + str(site_by_overhang["TGA"]["overhang_5to3"])
    if payload != supplemental:
        raise ValueError(
            "The linker-mNeonGreen-stop payload reconstructed from the uploaded backbone "
            "does not match the sequence supplied in Bollen supplementary S2."
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


def load_c_terminal_payload() -> str:
    """Return the payload extracted directly from the uploaded .dna backbone."""
    return fixed_backbone_analysis().payload_sequence


def backbone_metadata() -> dict[str, Any]:
    analysis = fixed_backbone_analysis()
    document = analysis.document
    return {
        "name": BACKBONE_NAME,
        "addgene_id": BACKBONE_ADDGENE_ID,
        "snapgene_file": document.path.name,
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
        "payload_stop_overhang_5to3": SAPI_OVERHANGS["payload_to_dha"],
        "payload_matches_bollen_s2": analysis.payload_matches_supplement,
    }


def payload_metadata() -> dict[str, str | int | bool]:
    payload = load_c_terminal_payload()
    linker_coding = payload[:21]
    tag_coding = payload[21:-3]
    payload_coding = payload[:-3]
    analysis = fixed_backbone_analysis()
    return {
        "name": "GGGGSAS linker + mNeonGreen + stop",
        "backbone_name": BACKBONE_NAME,
        "addgene_id": BACKBONE_ADDGENE_ID,
        "payload_sequence_5to3": payload,
        "payload_length_nt": len(payload),
        "linker_coding_sequence": linker_coding,
        "linker_peptide": translate(linker_coding),
        "tag_coding_sequence": tag_coding,
        "tag_length_nt": len(tag_coding),
        "tag_peptide": translate(tag_coding),
        "tag_length_aa": len(translate(tag_coding)),
        "payload_coding_sequence": payload_coding,
        "payload_coding_length_nt": len(payload_coding),
        "payload_peptide": translate(payload_coding),
        "payload_peptide_length_aa": len(translate(payload_coding)),
        "stop_codon": payload[-3:],
        "source": (
            "Extracted from the uploaded Addgene #169227 SnapGene file and verified "
            "base-for-base against Bollen et al. supplementary S2"
        ),
        "matches_bollen_s2": analysis.payload_matches_supplement,
    }


def _shifted_backbone_features(
    analysis: FixedBackboneAnalysis, *, replacement_length: int
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
                "note": "Retained from uploaded Addgene #169227 SnapGene annotation",
            }
        )
    return retained


def assemble_c_terminal_plasmid(
    *, target_with_pam: str, uha: str, dha: str
) -> dict[str, Any]:
    """Simulate the SapI Golden Gate product in uploaded Addgene #169227."""
    target = clean_dna(target_with_pam)
    uha = clean_dna(uha)
    dha = clean_dna(dha)
    if len(target) != 23:
        raise ValueError("SpCas9 target-with-PAM must be 23 nt")
    if any(motif in uha or motif in dha for motif in SAPI_RECOGNITION_MOTIFS):
        raise ValueError("Final homology arms must be SapI-free before plasmid assembly")

    analysis = fixed_backbone_analysis()
    payload = analysis.payload_sequence
    donor_insert = target + uha + payload + dha + target
    original = analysis.document.sequence

    # Simulate the top-strand Golden Gate product explicitly. The TAC cohesive
    # end is retained once before the 5-prime donor target, while AAT is already
    # the first three bases of the retained right-hand vector fragment.
    assembled_replacement = SAPI_OVERHANGS["vector_to_uha"] + donor_insert
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

    donor_start0 = analysis.vector_left_cut0 + len(SAPI_OVERHANGS["vector_to_uha"])
    target_5_start0 = donor_start0
    uha_start0 = target_5_start0 + len(target)
    payload_start0 = uha_start0 + len(uha)
    linker_end0 = payload_start0 + 21
    tag_start0 = linker_end0
    payload_end0 = payload_start0 + len(payload)
    dha_start0 = payload_end0
    target_3_start0 = dha_start0 + len(dha)
    donor_end0 = target_3_start0 + len(target)

    # The junction overhangs are retained once in the ligated top-strand sequence.
    junctions = {
        "vector_to_uha": {
            "overhang_5to3": "TAC",
            "start0": donor_start0 - 3,
            "end0": donor_start0,
            "observed": assembled[donor_start0 - 3:donor_start0],
            "window_5to3": assembled[donor_start0 - 12:donor_start0 + 12],
        },
        "uha_to_payload": {
            "overhang_5to3": "GGC",
            "start0": payload_start0,
            "end0": payload_start0 + 3,
            "observed": assembled[payload_start0:payload_start0 + 3],
            "window_5to3": assembled[payload_start0 - 12:payload_start0 + 12],
        },
        "payload_to_dha": {
            "overhang_5to3": "TGA",
            "start0": payload_end0 - 3,
            "end0": payload_end0,
            "observed": assembled[payload_end0 - 3:payload_end0],
            "window_5to3": assembled[payload_end0 - 12:payload_end0 + 12],
        },
        "dha_to_vector": {
            "overhang_5to3": "AAT",
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
    features = _shifted_backbone_features(analysis, replacement_length=len(assembled_replacement))
    features.extend(
        [
            {
                "type": "misc_feature",
                "label": "assembled donor cassette",
                "start0": donor_start0,
                "end0": donor_end0,
                "strand": 1,
                "note": "target-UHA-linker-mNeonGreen-stop-DHA-target",
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
                "label": "GGGGSAS linker",
                "start0": payload_start0,
                "end0": linker_end0,
                "strand": 1,
                "note": "First GGC codon also forms the UHA-to-payload SapI overhang",
            },
            {
                "type": "CDS",
                "label": "mNeonGreen + stop",
                "start0": tag_start0,
                "end0": payload_end0,
                "strand": 1,
                "note": "Fixed Addgene #169227 mNeonGreen coding sequence including TGA stop",
            },
            {
                "type": "misc_feature",
                "label": "3-prime homology arm (DHA), final",
                "start0": dha_start0,
                "end0": target_3_start0,
                "strand": 1,
                "note": "Gene-oriented final arm; endogenous stop omitted",
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
        "uploaded_backbone": backbone_metadata(),
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
            "linker_end0": linker_end0,
            "tag_start0": tag_start0,
            "payload_end0": payload_end0,
            "dha_start0": dha_start0,
            "target_3_start0": target_3_start0,
        },
    }


def c_terminal_synthesis_fragments(
    *, target_with_pam: str, uha: str, dha: str
) -> dict[str, Any]:
    target = clean_dna(target_with_pam)
    uha = clean_dna(uha)
    dha = clean_dna(dha)
    if len(target) != 23:
        raise ValueError("SpCas9 target-with-PAM must be 23 nt")
    uha_fragment = CTERM_UHA_PREFIX + target + uha + CTERM_UHA_SUFFIX
    dha_fragment = CTERM_DHA_PREFIX + dha + target + CTERM_DHA_SUFFIX
    uha_sapi_sites = find_sapi_sites(uha_fragment)
    dha_sapi_sites = find_sapi_sites(dha_fragment)
    if [site["overhang_5to3"] for site in uha_sapi_sites] != ["TAC", "GGC"]:
        raise ValueError("The UHA synthesis fragment does not generate TAC/GGC SapI overhangs")
    if [site["overhang_5to3"] for site in dha_sapi_sites] != ["TGA", "AAT"]:
        raise ValueError("The DHA synthesis fragment does not generate TGA/AAT SapI overhangs")
    assembly = assemble_c_terminal_plasmid(
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
        "expected_sapi_overhangs": dict(SAPI_OVERHANGS),
        **assembly,
    }
