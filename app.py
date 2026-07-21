from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import tempfile
from typing import Any

import pandas as pd
import streamlit as st

from hdr_designer.design import DesignError, design_online, design_tubb5_fixture
from hdr_designer.backbones import backbone_for_terminus, infer_custom_backbone_definition
from hdr_designer.ensembl import EnsemblError, SPECIES
from hdr_designer.exports import (
    arms_fasta,
    assembled_plasmid_genbank,
    design_json,
    design_report,
    guide_rows,
    guides_csv,
    sapi_qc_rows,
)
from hdr_designer.models import DesignResult, HomologyArm
from hdr_designer.snapgene import SnapGeneError

APP_VERSION = "0.5.0"


def _download_buttons(result: DesignResult) -> None:
    stem = f"{result.gene_symbol}_{result.terminus.lower().replace('-', '_').replace(' ', '_')}"
    cols = st.columns(5)
    cols[0].download_button(
        "Full report (.txt)", design_report(result), f"{stem}_design_report.txt", "text/plain",
        width="stretch",
    )
    cols[1].download_button(
        "Guides (.csv)", guides_csv(result), f"{stem}_guides.csv", "text/csv",
        width="stretch",
    )
    cols[2].download_button(
        "Sequences (.fasta)", arms_fasta(result), f"{stem}_sequences.fasta", "text/plain",
        width="stretch",
    )
    cols[3].download_button(
        "Assembled plasmid (.gb)",
        assembled_plasmid_genbank(result) if result.sequence_complete else "",
        f"{stem}_assembled_plasmid.gb",
        "text/plain",
        width="stretch",
        disabled=not result.sequence_complete,
    )
    cols[4].download_button(
        "Full data (.json)", design_json(result), f"{stem}_design.json", "application/json",
        width="stretch",
    )


def _show_arm(arm: HomologyArm) -> None:
    st.dataframe(
        pd.DataFrame(
            [{
                "Arm": arm.name,
                "Reference interval": arm.genomic_interval_1based,
                "Length (bp)": arm.length,
                "GC (%)": arm.gc_percent,
                "Raw SapI sites": len(arm.sapi_sites),
                "Final SapI sites": len(arm.final_sapi_sites),
                "Verified synonymous changes": len(arm.mutations),
            }]
        ),
        hide_index=True,
        width="stretch",
    )
    if arm.mutations:
        st.markdown("**Verified synonymous changes in final arm**")
        st.dataframe(
            pd.DataFrame([asdict(mutation) for mutation in arm.mutations]),
            hide_index=True,
            width="stretch",
        )
    with st.expander("Reference gene-oriented sequence"):
        st.code(arm.gene_oriented_sequence, language=None)
    if arm.final_gene_oriented_sequence != arm.gene_oriented_sequence:
        with st.expander("FINAL gene-oriented sequence", expanded=True):
            st.code(arm.final_gene_oriented_sequence, language=None)
    with st.expander("Chromosome-forward reference sequence"):
        st.code(arm.chromosome_forward_sequence, language=None)


def _show_sapi_quality_control(result: DesignResult) -> None:
    st.markdown("### SapI arm quality control")
    arms = (result.five_prime_arm, result.three_prime_arm)
    rows = sapi_qc_rows(result)
    found = sum(len(arm.sapi_sites) for arm in arms)
    remaining = sum(len(arm.final_sapi_sites) for arm in arms)
    resolved = sum(row["Status"] == "Resolved" for row in rows)
    cols = st.columns(3)
    cols[0].metric("Arm SapI sites found", found)
    cols[1].metric("Arm SapI sites resolved", resolved)
    cols[2].metric("Arm SapI sites remaining", remaining)

    arm_rows = []
    for arm in arms:
        arm_resolved = sum(
            row["Arm"] == arm.name and row["Status"] == "Resolved"
            for row in rows
        )
        arm_rows.append(
            {
                "Arm": arm.name,
                "Sites found": len(arm.sapi_sites),
                "Sites resolved": arm_resolved,
                "Sites remaining": len(arm.final_sapi_sites),
            }
        )
    st.dataframe(pd.DataFrame(arm_rows), hide_index=True, width="stretch")

    if found == 0 and remaining == 0:
        st.success("No internal SapI recognition site was found in either homology arm.")
    elif remaining == 0:
        st.success(
            f"All {found} internal SapI site(s) were removed from the final arm sequences."
        )
    else:
        st.warning(
            f"{remaining} SapI site(s) remain in the final arms; order-ready cloning "
            "fragments are withheld pending review."
        )
    if rows:
        st.markdown("#### Site-by-site resolution")
        st.dataframe(pd.DataFrame(rows), hide_index=True, width="stretch")


def _show_result(result: DesignResult) -> None:
    st.divider()
    if result.sequence_complete:
        st.success(result.status)
    else:
        st.warning(result.status)

    st.subheader(f"{result.gene_symbol}: {result.terminus} design")
    cols = st.columns(6)
    cols[0].metric("Assembly", result.assembly)
    cols[1].metric("Transcript", result.transcript_id)
    cols[2].metric("Guides", len(result.guides))
    cols[3].metric("Arms", f"{result.homology_arm_length} bp")
    cols[4].metric("Native protein", f"{result.protein_length_aa} aa")
    cols[5].metric("Fusion", f"{result.fusion_protein_length_aa or 0} aa")
    st.caption(
        f"chr{result.chromosome}, gene strand {result.gene_strand}; insertion boundary "
        f"between 1-based bases {result.insertion_boundary0:,} and {result.insertion_boundary0 + 1:,}."
    )

    st.markdown("### Validation gates")
    st.dataframe(pd.DataFrame(result.validations), hide_index=True, width="stretch")

    if result.guides:
        guide = result.top_guide
        st.markdown("### Selected guide")
        cols = st.columns(6)
        cols[0].metric("Spacer", guide.spacer)
        cols[1].metric("PAM", guide.pam)
        cols[2].metric("PAM strand", guide.chromosome_strand)
        cols[3].metric("Nick distance", f"{guide.distance_to_insertion} bp")
        cols[4].metric("Final PAM disrupted", "yes" if guide.final_pam_destroyed else "no")
        cols[5].metric("Retained segment", f"{guide.final_longest_retained_segment} nt")
        st.code(f"5'-{guide.spacer}-{guide.pam}-3'", language=None)
        st.write(guide.rationale)
        if guide.blocking_mutation_note:
            st.info(guide.blocking_mutation_note)
        st.caption(result.guide_scoring_note)
        with st.expander("All ranked candidates"):
            st.dataframe(pd.DataFrame(guide_rows(result)), hide_index=True, width="stretch")

    _show_sapi_quality_control(result)

    st.markdown("### Homology arms")
    tab_uha, tab_dha = st.tabs(["5-prime arm (UHA)", "3-prime arm (DHA)"])
    with tab_uha:
        _show_arm(result.five_prime_arm)
    with tab_dha:
        _show_arm(result.three_prime_arm)

    if result.donor_payload:
        st.markdown("### Fixed donor payload")
        payload_cols = st.columns(4)
        payload_cols[0].metric("Linker", str(result.donor_payload.get("linker_peptide", "")))
        payload_cols[1].metric(
            str(result.donor_payload.get("tag_name", "Tag")),
            f"{result.donor_payload.get('tag_length_aa', 0)} aa",
        )
        payload_cols[2].metric("Payload coding", f"{result.donor_payload.get('payload_coding_length_nt', 0)} nt")
        payload_cols[3].metric("Stop", str(result.donor_payload.get("stop_codon", "")))
        with st.expander("Payload sequence"):
            st.code(str(result.donor_payload.get("payload_sequence_5to3", "")), language=None)

    if result.cloning_fragments:
        st.markdown("### SapI / Golden Gate fragments")
        overhangs = result.cloning_fragments.get("expected_sapi_overhangs")
        if overhangs:
            st.write("Expected overhangs:", overhangs)
        fragment_sites = result.cloning_fragments.get("synthesis_fragment_sapi_sites")
        if fragment_sites:
            st.write(
                "Verified fragment digests:",
                {
                    name: [site.get("overhang_5to3") for site in sites]
                    for name, sites in fragment_sites.items()
                },
            )
        for label, key in (
            ("FINAL UHA synthesis fragment", "uha_synthesis_fragment_5to3"),
            ("FINAL DHA synthesis fragment", "dha_synthesis_fragment_5to3"),
            ("Assembled donor insert simulation", "assembled_donor_insert_5to3"),
            ("UHA preview", "uha_synthesis_fragment_preview_5to3"),
            ("DHA preview", "dha_synthesis_fragment_preview_5to3"),
        ):
            seq = result.cloning_fragments.get(key)
            if seq:
                with st.expander(f"{label} ({len(str(seq))} nt)", expanded="FINAL" in label):
                    st.code(str(seq), language=None)

        backbone = result.cloning_fragments.get("uploaded_backbone")
        plasmid = result.cloning_fragments.get("assembled_plasmid_5to3")
        if backbone and plasmid:
            st.markdown("#### Uploaded backbone and full in-silico assembly")
            cols = st.columns(5)
            cols[0].metric("Input backbone", f"{backbone.get('length_nt')} bp")
            cols[1].metric("Input SapI sites", str(backbone.get("sapi_site_count")))
            cols[2].metric("Final plasmid", f"{result.cloning_fragments.get('assembled_plasmid_length_nt')} bp")
            cols[3].metric("Final SapI sites", str(result.cloning_fragments.get("assembled_plasmid_sapi_site_count")))
            cols[4].metric("Payload verified", "yes" if backbone.get("payload_sequence_verified") else "no")
            st.caption(
                f"Parsed {backbone.get('snapgene_file')} | SHA-256 {backbone.get('snapgene_sha256')}"
            )
            junction_rows = []
            for name, details in result.cloning_fragments.get("golden_gate_junctions", {}).items():
                junction_rows.append(
                    {
                        "Junction": name,
                        "Expected overhang": details.get("overhang_5to3"),
                        "Observed": details.get("observed"),
                        "Sequence window (5'->3')": details.get("window_5to3"),
                    }
                )
            if junction_rows:
                st.dataframe(pd.DataFrame(junction_rows), hide_index=True, width="stretch")
            with st.expander(f"Full assembled circular plasmid sequence ({len(str(plasmid))} bp)"):
                st.code(str(plasmid), language=None)

    if result.primer_tail_templates:
        st.markdown("### PCR-primer 5-prime tail templates")
        st.dataframe(
            pd.DataFrame(
                [{"Template": key, "Sequence / note": value}
                 for key, value in result.primer_tail_templates.items()]
            ),
            hide_index=True,
            width="stretch",
        )

    st.markdown("### Export")
    _download_buttons(result)

    with st.expander("Warnings and provenance"):
        st.markdown("**Warnings**")
        for warning in result.warnings:
            st.write(f"- {warning}")
        st.markdown("**Provenance**")
        for item in result.provenance:
            st.write(f"- {item}")
        preview: dict[str, Any] = asdict(result)
        for arm_key in ("five_prime_arm", "three_prime_arm"):
            preview[arm_key]["gene_oriented_sequence"] = "[sequence in downloads]"
            preview[arm_key]["chromosome_forward_sequence"] = "[sequence in downloads]"
            preview[arm_key]["corrected_gene_oriented_sequence"] = "[sequence in downloads]"
        for key in (
            "uha_synthesis_fragment_5to3",
            "dha_synthesis_fragment_5to3",
            "assembled_donor_insert_5to3",
            "assembled_plasmid_5to3",
        ):
            if key in preview.get("cloning_fragments", {}):
                preview["cloning_fragments"][key] = "[sequence in downloads]"
        if preview.get("donor_payload"):
            preview["donor_payload"]["payload_sequence_5to3"] = "[sequence in downloads]"
            preview["donor_payload"]["tag_coding_sequence"] = "[sequence in downloads]"
            preview["donor_payload"]["payload_coding_sequence"] = "[sequence in downloads]"
        st.json(preview)


def main() -> None:
    st.set_page_config(page_title="HDR Tag Designer", page_icon="DNA", layout="wide")
    st.title("HDR Tag Designer")
    st.caption(
        f"Bollen-style ITPN gene-tagging prototype using SpCas9 D10A. Version {APP_VERSION} "
        "adds complete N- and C-terminal transfer-vector assembly."
    )
    st.info(
        "Species is the first design choice. The bundled Tubb5 test is reproducible offline; "
        "other human or mouse genes require live Ensembl access. No off-target analysis is run."
    )

    with st.sidebar:
        st.header("Design inputs")
        species_key = st.selectbox(
            "1. Species",
            ["mouse", "human"],
            format_func=lambda key: f"{SPECIES[key].label} - {SPECIES[key].assembly}",
            index=0,
        )
        terminus = st.selectbox("2. Tagging terminus", ["C-terminal", "N-terminal"], index=0)
        if species_key == "mouse":
            if not terminus.startswith("C"):
                st.session_state["use_tubb5_fixture"] = False
            use_fixture = st.checkbox(
                "Use bundled Tubb5 validation fixture",
                value=terminus.startswith("C"),
                key="use_tubb5_fixture",
                disabled=not terminus.startswith("C"),
                help=(
                    "Offline C-terminal Tubb5-201 test with 600-bp arms. "
                    "The bundled transcript is ENSMUST00000001566.10."
                ),
            ) and terminus.startswith("C")
        else:
            use_fixture = False
            st.caption("The bundled offline fixture is mouse-only; human designs use live Ensembl.")

        gene = st.text_input(
            "3. Gene symbol or Ensembl gene ID",
            value="Tubb5" if species_key == "mouse" else "TUBB5",
            key=f"gene_{species_key}",
            disabled=use_fixture,
        )
        transcript_id = st.text_input(
            "4. Transcript ID (optional)",
            value="ENSMUST00000001566.10" if use_fixture else "",
            key=f"transcript_{species_key}_{'fixture' if use_fixture else 'live'}",
            disabled=use_fixture,
            help=(
                "Bundled fixture version; this is the current canonical stable transcript listing checked for the build."
                if use_fixture
                else "Blank uses the live Ensembl canonical transcript."
            ),
        )
        selected_backbone = backbone_for_terminus(terminus)
        backbone_mode = st.radio(
            "5. Donor backbone",
            ["Verified built-in", "Upload custom .dna"],
            horizontal=True,
        )
        custom_backbone_upload = None
        if backbone_mode == "Verified built-in":
            st.text_input(
                "Selected backbone",
                value=f"{selected_backbone.name} (Addgene #{selected_backbone.addgene_id})",
                disabled=True,
            )
        else:
            custom_backbone_upload = st.file_uploader(
                "Custom SnapGene backbone",
                type=["dna"],
                help=(
                    "Must be circular and retain exactly four SapI sites in the Bollen "
                    "N- or C-terminal overhang order plus the GGGGSAS linker."
                ),
            )
        arm_length = st.number_input(
            "6. Homology arm length (bp)", 100, 2000, 600, 50
        )
        guide_window = st.number_input(
            "7. Guide nick search radius (bp)", 10, 200, 50, 10
        )
        run = st.button("Design locus", type="primary", width="stretch")

    st.markdown("#### Scope of this version")
    st.write(
        "ITPN / SpCas9 D10A; mouse GRCm39 or human GRCh38; reference sequence only; "
        "verified Addgene #169226 N-terminal and #169227 C-terminal mNeonGreen backbones, "
        "plus structurally validated custom SnapGene backbones; "
        "automatic synonymous coding mutations with noncoding cases held for review."
    )

    if not run:
        st.write("Use the default settings and click **Design locus** to run the Tubb5 test.")
        return

    temporary_backbone_path: Path | None = None
    try:
        with st.spinner("Building and validating design..."):
            design_backbone = selected_backbone
            if backbone_mode == "Upload custom .dna":
                if custom_backbone_upload is None:
                    raise DesignError("Upload a custom SnapGene .dna backbone before designing.")
                with tempfile.NamedTemporaryFile(
                    prefix="hdr_tag_custom_", suffix=".dna", delete=False
                ) as handle:
                    handle.write(custom_backbone_upload.getvalue())
                    temporary_backbone_path = Path(handle.name)
                design_backbone = infer_custom_backbone_definition(temporary_backbone_path)
                if design_backbone.terminus != terminus:
                    raise DesignError(
                        f"The uploaded SapI architecture is {design_backbone.terminus}, "
                        f"but {terminus} tagging was selected."
                    )
            if use_fixture:
                if backbone_mode != "Verified built-in":
                    raise DesignError("The bundled Tubb5 fixture uses the verified built-in backbone.")
                if species_key != "mouse" or gene.strip().casefold() != "tubb5":
                    raise DesignError("The bundled fixture is fixed to mouse Tubb5.")
                if not terminus.startswith("C") or int(arm_length) != 600:
                    raise DesignError("The bundled fixture is fixed to C-terminal tagging with 600-bp arms.")
                result = design_tubb5_fixture(
                    arm_length=int(arm_length), guide_window=int(guide_window)
                )
            else:
                result = design_online(
                    species_key=species_key,
                    gene=gene,
                    transcript_id=transcript_id or None,
                    terminus=terminus,
                    arm_length=int(arm_length),
                    guide_window=int(guide_window),
                    backbone_definition=design_backbone,
                )
    except (DesignError, EnsemblError, SnapGeneError, ValueError) as exc:
        st.error(str(exc))
        st.stop()
    finally:
        if temporary_backbone_path is not None and temporary_backbone_path.exists():
            temporary_backbone_path.unlink()

    st.session_state["latest_design"] = result
    _show_result(result)


if __name__ == "__main__":
    main()
