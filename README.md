# HDR Tag Designer

A local Streamlit prototype for Bollen-style **in-trans paired nicking (ITPN)** gene tagging with **SpCas9 D10A**. The first input is species: mouse **GRCm39** or human **GRCh38**.

Version 0.5.1 supports both finalized Bollen donor architectures:

- N-terminal `mNeonGreen-GGGGSAS` using Addgene **#169226**
- C-terminal `GGGGSAS-mNeonGreen-stop` using Addgene **#169227**
- uploaded custom SnapGene backbones that retain one of those two four-SapI/linker architectures
- 600-bp arms by default
- no off-target analysis
- reference sequence only

Custom files are accepted only after exact structural checks; arbitrary adapters, linkers, enzymes, and payload-only FASTA uploads remain deferred.

## What it does

1. Resolves a gene and transcript through public Ensembl REST in live mode.
2. Defines an N- or C-terminal insertion boundary from the selected transcript.
3. Finds all SpCas9-NGG candidates with a nominal nick within the chosen window.
4. Ranks guides by the Bollen priorities available locally: distance first, target disruption second, then basic GC/poly-T properties. A nearer guide is not demoted merely because it needs a synonymous blocking mutation.
5. Generates gene-oriented homology arms, detects every internal SapI site, and automatically removes it using either a verified synonymous coding change or a guarded single-base non-coding change. Non-coding candidates exclude coding bases, the first/last three bases of every exon, and six intronic bases on either side of each splice boundary; they must not create another SapI site or a longer problematic homopolymer. The automatic ranking prefers a base outside the mature transcript, lower local homopolymer burden, a central motif position, then a transition.
6. Re-evaluates the selected target after all donor edits. If needed, it searches synonymous changes that destroy the PAM first, then PAM-proximal seed changes that satisfy the 14-nt retained-segment cutoff. Every released change must preserve the complete CDS translation and avoid new SapI sites.
7. Generates the exact Bollen supplementary N- or C-terminal UHA/DHA synthesis fragments and PCR-primer tail templates.
8. Parses Addgene #169226 or #169227, verifies all four SapI sites, and reconstructs the payload between the inner cuts. N-terminal order is `TAC -> GTG -> AGC -> AAT`; C-terminal order is `TAC -> GGC -> TGA -> AAT`.
9. Reconstructs the full circular Golden Gate product and verifies all four ligation junctions and the absence of residual SapI sites.
10. Exports TXT, CSV, FASTA, JSON, and an annotated GenBank file that can be opened in SnapGene and similar sequence editors.
11. Presents a SapI quality-control panel with total and per-arm counts plus a site-by-site record of coordinates, resolution status, nucleotide/codon changes, protein consequence, and selection reason.
12. Optionally classifies an uploaded circular SnapGene backbone as N- or C-terminal from its SapI overhang order, validates its GGGGSAS linker and reading frame, then runs the same complete assembly and export path.

The prototype does **not** calculate a validated on-target activity score, perform off-target searches, inspect sample-specific variants, or design locus-specific PCR annealing regions. Non-coding guide-blocking changes remain withheld, while non-coding SapI sites are now domesticated automatically when an eligible one-base solution passes the sequence gates. Sites confined to coding bases without a synonymous solution, protected splice-edge bases, or otherwise unsafe candidates remain blocked.

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

The fixed payload is 729 bp (`GGGGSAS` linker + 235-aa mNeonGreen + `TGA`). The predicted fusion is 686 aa. The donor insert is 1,975 bp and the simulated circular plasmid is 3,950 bp.

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

For a custom `.dna` upload, the app requires a circular record, exactly four SapI sites in one supported order, a frame-compatible payload, the expected GGGGSAS junction linker, and no internal in-frame stop. A cassette crossing the file's sequence origin must first be rotated so the `TAC` cut precedes the payload.

## Important limitations

This is a computational design aid, not an experimentally validated clinical or diagnostic tool. Independently verify transcript choice, all junctions, the current plasmid sequence, sample genotype, guide activity, and biological suitability before ordering or experimentation. C-terminal tagging of Tubb5 may perturb its functionally important tubulin tail; the bundled result is computational only.

## Method and sequence sources

- Bollen et al. (2022), paper and supplementary S1-S3: guide priorities, 14-nt retargeting cutoff, 600-bp arms, SapI adapters, and mNeonGreen payload.
- Uploaded Addgene plasmids #169226 and #169227: fixed circular backbone sequences and annotations.
- `data/bollen_supplementary_s1.docx`: N- and C-terminal fragment/primer architecture supplied from the original paper.
- Ensembl REST: live human/mouse transcript and reference-sequence retrieval.
- UCSC/GENCODE mm39: bundled Tubb5 reference fixture.
