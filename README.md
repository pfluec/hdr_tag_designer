# HDR Tag Designer

A local Streamlit prototype for Bollen-style **in-trans paired nicking (ITPN)** gene tagging with **SpCas9 D10A**. The first input is species: mouse **GRCm39** or human **GRCh38**.

Version 0.7.0 supports both finalized Bollen donor architectures:

- N-terminal `mNeonGreen-GGGGSAS` using Addgene **#169226**
- C-terminal `GGGGSAS-mNeonGreen-stop` using Addgene **#169227**
- uploaded custom SnapGene backbones that retain one of those two four-SapI architectures, including multi-ORF and non-frame payload cassettes
- 600-bp arms by default
- no off-target analysis
- reference sequence only

Custom files are accepted only after exact SapI structural checks. Complex cassettes are assembled without asserting a single fusion translation, and non-frame payloads receive an explicit warning. Arbitrary adapters, enzymes, and payload-only FASTA uploads remain deferred.

## What it does

1. Resolves a gene and transcript through public Ensembl REST in live mode.
2. Defines an N- or C-terminal insertion boundary from the selected transcript.
3. Finds all SpCas9-NGG candidates with a nominal nick within the chosen window.
4. Ranks guides by the Bollen priorities available locally: distance first, target disruption second, then basic GC/poly-T properties. A nearer guide is not demoted merely because it needs a synonymous blocking mutation.
5. Generates gene-oriented homology arms, detects every internal SapI site, and automatically removes it using either a verified synonymous coding change or a guarded single-base non-coding change. Non-coding candidates exclude coding bases, the first/last three bases of every exon, and six intronic bases on either side of each splice boundary; they must not create another SapI site or a longer problematic homopolymer. The automatic ranking prefers a base outside the mature transcript, lower local homopolymer burden, a central motif position, then a transition.
6. Reconstructs the complete edited allele (`final UHA + payload + final DHA`) and scans both strands for contiguous 20-nt protospacers immediately adjacent to NGG. A surviving or junction-recreated on-target candidate requires blocking only when its contiguous PAM-proximal match to the original spacer exceeds the 14-nt cutoff. If needed, the tool searches synonymous PAM changes first and then synonymous seed changes; every released change must preserve the complete CDS translation, avoid new SapI sites, and eliminate every recuttable candidate from the final edited locus. If the highest-ranked guide cannot pass this gate, the next ranked sequence-complete guide is tried.
7. Generates the exact Bollen supplementary N- or C-terminal UHA/DHA fragments and complete homology-arm cloning primers, combining the fixed SapI/Golden Gate tails with Primer3-evaluated genomic annealing regions.
8. Parses Addgene #169226 or #169227, verifies all four SapI sites, and reconstructs the payload between the inner cuts. N-terminal order is `TAC -> GTG -> AGC -> AAT`; C-terminal order is `TAC -> GGC -> TGA -> AAT`.
9. Reconstructs the full circular Golden Gate product and verifies all four ligation junctions and the absence of residual SapI sites.
10. Shows the reference guide target/PAM, the target after automatic point mutations, and the actual donor-allele sequence across the edited region, including an insertion marker and the full inserted context.
11. Designs WT-locus, 5-prime junction, and 3-prime junction genotyping PCR assays with Primer3 thermodynamics. Junction genomic primers bind outside the homology arms, payload primers bind at least 150 bp from the tested junction, and external primers are checked against the assembled donor plasmid.
12. Exports one deterministic ordering ZIP containing a Twist sequence CSV, the two selected-guide cloning oligos as CSV, four unique genotyping primers with full assay metadata as CSV, and annotated GenBank records for the assembled plasmid, WT locus, and edited locus. The older report/FASTA/JSON generators remain available internally but are no longer separate UI downloads.
13. Presents a SapI quality-control panel with total and per-arm counts plus a site-by-site record of coordinates, resolution status, nucleotide/codon changes, protein consequence, and selection reason.
14. Optionally classifies an uploaded circular SnapGene backbone as N- or C-terminal from its SapI overhang order, then runs the same complete assembly and export path. A conventional GGGGSAS single-ORF payload receives a fusion prediction; other payload organizations are retained as complex cassettes with a warning.
15. Builds separate annotated linear WT and edited locus records, each extending 300 bp beyond both homology arms. The records annotate reference/final arms, payload or native junction sequence, and applicable genotyping-primer binding sites and are included in the ordering ZIP as GenBank.
16. Applies an ordering preflight to the final, mutation-containing homology arms. If an arm contains a homopolymer longer than 14 nt, the tool automatically moves the distal arm boundary into the middle of the relevant run while preserving the insertion-proximal homology. The retained partial run sits at the guide/SapI-facing arm edge and is at most 14 nt. The requested/final arm lengths and exact adjustment are reported in the UI, CSV, JSON, reports, and GenBank annotations. Twist portal screening is still required for all submitted sequences.

The prototype does **not** calculate a validated on-target activity score, perform off-target searches, inspect sample-specific variants, or perform a genome-wide primer-uniqueness search. Genotyping primers are checked against the donor plasmid, but must still be confirmed with Primer-BLAST before ordering. Non-coding guide-blocking changes remain withheld, while non-coding SapI sites are now domesticated automatically when an eligible one-base solution passes the sequence gates. Sites confined to coding bases without a synonymous solution, protected splice-edge bases, or otherwise unsafe candidates remain blocked.

There is not yet an interactive mutation editor. Accordingly, cases the automatic engine cannot resolve are labeled `DESIGN BLOCKED - NO SEQUENCE-COMPLETE OUTPUT` rather than suggesting that an in-app review action is available. A future mutation-choice interface can expose alternative eligible SapI substitutions.

## Run locally

```zsh
conda env create -f environment.yml
conda activate hdr-tag-designer
python -m pip install -e .
streamlit run app.py
```

After the environment exists, `./run.zsh` activates it, refreshes the runtime dependencies, and starts Streamlit.

## Run the tests and regenerate the bundled result

```zsh
./run_tests.zsh
```

Generated Tubb5 files are written to `outputs/tubb5_test/`.

## Bundled Tubb5 validation

The offline test uses mouse **Tubb5-201 / ENSMUST00000001566.10** on GRCm39/mm39. It validates the 444-aa coding endpoint and the terminal reference-genome interval used for the design. A live Ensembl check on 2026-07-21 resolved the same stable transcript as version `.11`; the fixture version is retained for reproducibility.

The selected C-terminal guide is:

```text
5'-GAGGCAGAAGAGGAGGCCTA-AGG-3'
```

Its nominal D10A nick lies 1 bp from the insertion boundary. Stop-codon replacement deletes part of the guide's `AGG` PAM, so the complete target is absent from the edited locus and no additional guide-blocking mutation is required. One synonymous `GAG -> GAA` change removes an internal SapI site from the UHA.

The Tubb5 fixture keeps its sequence-locked regression correction. Live mouse and human designs use the generic synonymous mutation search and report each genomic/arm coordinate, codon consequence, PAM change, retained-target change, and selection reason in TXT/JSON/GenBank outputs.

The fixed payload is 729 bp (`GGGGSAS` linker + 235-aa mNeonGreen + `TGA`). The predicted fusion is 686 aa. After the homopolymer-aware DHA adjustment, the donor insert is 1,643 bp and the simulated circular plasmid is 3,618 bp.

The bundled result also locks three genotyping assays: a 1,014-bp WT-locus product (1,740 bp from the edited allele), an 845-bp 5-prime junction product, and a 528-bp 3-prime junction product. Both external locus primers are outside the final homology arms and absent from the assembled donor plasmid; both payload primers are more than 150 bp from their tested junction.

The requested 600-bp Tubb5 3-prime homology arm contains a 26-nt poly-T run at original arm bases 256-281. The tool retains the insertion-proximal 268 bp and terminates the DHA after the first 13 T bases, immediately before the external guide/SapI architecture. The final UHA/DHA lengths are therefore 600/268 bp, the order-ready ZIP passes local homopolymer QC, and no genomic base is mutated to resolve the synthesis constraint.

The four exported genotyping oligos are named `Tubb5_wt_5_fwd`, `Tubb5_wt_3_rev`, `Tubb5_mut_5_rev`, and `Tubb5_mut_3_fwd`. The external WT primers are reused exactly in their corresponding 5-prime and 3-prime junction assays instead of being emitted twice.

## Uploaded-backbone verification

The bundled `data/addgene_169227.dna` file is parsed directly. The expected input properties are:

- length: 2,768 bp
- topology: circular
- SHA-256: `b5ccaa5a257b71a1f2bac05ab15785f07098a46ed805c2862b5beeade04046b1`
- SapI recognition sites: four
- top-strand overhangs in assembly order: `TAC`, `GGC`, `TGA`, `AAT`
- extracted linker-mNeonGreen-stop payload: 729 bp and identical to the Bollen S2 sequence

The assembled GenBank file includes the retained backbone annotations, the two donor target sites, both homology arms, linker, mNeonGreen CDS and stop, the four SapI junctions, and the single synonymous Tubb5 SapI-domestication change.

The bundled `data/addgene-169226.dna` N-terminal backbone is also verified directly:

- length/topology: 2,765 bp / circular
- SHA-256: `75d9d25b4dac8083c401ee5ac76b080a5f62e42d23357a81b4ac33e84a434177`
- SapI overhang order: `TAC`, `GTG`, `AGC`, `AAT`
- extracted payload: 726 bp (235-aa mNeonGreen followed by GGGGSAS; no stop codon)
- 600-bp-arm donor insert/final plasmid: 1,972 bp / 3,947 bp

For a custom `.dna` upload, the app requires a circular record and exactly four SapI sites in one supported order. A conventional frame-compatible payload with the expected GGGGSAS junction linker receives the single-fusion interpretation. Multi-ORF, non-coding, internal-stop, missing-linker, and non-frame cassettes are accepted and assembled, but are labelled as complex payloads and no single fusion-protein translation is asserted. A cassette crossing the file's sequence origin must first be rotated so the `TAC` cut precedes the payload.

## Important limitations

This is a computational design aid, not an experimentally validated clinical or diagnostic tool. Independently verify transcript choice, all junctions, the current plasmid sequence, sample genotype, guide activity, and biological suitability before ordering or experimentation. C-terminal tagging of Tubb5 may perturb its functionally important tubulin tail; the bundled result is computational only.

## Method and sequence sources

- Bollen et al. (2022), paper and supplementary S1-S3: guide priorities, 14-nt retargeting cutoff, 600-bp arms, SapI adapters, and mNeonGreen payload.
- Uploaded Addgene plasmids #169226 and #169227: fixed circular backbone sequences and annotations.
- `data/bollen_supplementary_s1.docx`: N- and C-terminal fragment/primer architecture supplied from the original paper.
- Ensembl REST: live human/mouse transcript and reference-sequence retrieval.
- UCSC/GENCODE mm39: bundled Tubb5 reference fixture.
- Primer3/primer3-py: PCR-primer melting temperatures and secondary-structure calculations; default 18-27-nt, 57-63 C, 35-65% GC rule set.
- Twist ordering platform: current CSV uploads use flexible name/sequence column mapping; the export therefore places `Name` and `Sequence` first and retains project QC metadata in additional columns. The local 14-nt homopolymer maximum is the project rule supplied by the user, while the Twist portal remains authoritative for complete synthesis screening.
