# HDR Tag Designer

A local Streamlit prototype for Bollen-style **in-trans paired nicking (ITPN)** gene tagging with **SpCas9 D10A**. The first input is species: mouse **GRCm39** or human **GRCh38**.

Version 0.4.0 deliberately supports one finalized donor architecture:

- C-terminal `GGGGSAS-mNeonGreen-stop`
- Bollen TVBB C-term-mNeongreen
- Addgene plasmid **#169227**
- the supplied Addgene `.dna` sequence as the fixed backbone reference
- 600-bp arms by default
- no off-target analysis
- reference sequence only

Custom donor backbones and custom tags remain deferred while the sequence-critical mutation workflow is validated.

## What it does

1. Resolves a gene and transcript through public Ensembl REST in live mode.
2. Defines an N- or C-terminal insertion boundary from the selected transcript.
3. Finds all SpCas9-NGG candidates with a nominal nick within the chosen window.
4. Ranks guides by the Bollen priorities available locally: distance first, target disruption second, then basic GC/poly-T properties. A nearer guide is not demoted merely because it needs a synonymous blocking mutation.
5. Generates gene-oriented homology arms, detects every internal SapI site, and automatically removes coding sites when a verified synonymous codon replacement is available.
6. Re-evaluates the selected target after all donor edits. If needed, it searches synonymous changes that destroy the PAM first, then PAM-proximal seed changes that satisfy the 14-nt retained-segment cutoff. Every released change must preserve the complete CDS translation and avoid new SapI sites.
7. Generates the exact Bollen S1/S3 C-terminal UHA and DHA synthesis fragments and PCR-primer tail templates.
8. Parses the uploaded Addgene #169227 SnapGene file, verifies its four SapI sites and `TAC -> GGC -> TGA -> AAT` overhang order, and checks that its 729-bp payload matches Bollen supplementary S2.
9. Reconstructs the full circular Golden Gate product and verifies all four ligation junctions and the absence of residual SapI sites.
10. Exports TXT, CSV, FASTA, JSON, and an annotated GenBank file that can be opened in SnapGene and similar sequence editors.
11. Presents a SapI quality-control panel with total and per-arm counts plus a site-by-site record of coordinates, resolution status, nucleotide/codon changes, protein consequence, and selection reason.

The prototype does **not** calculate a validated on-target activity score, perform off-target searches, inspect sample-specific variants, or design locus-specific PCR annealing regions. Noncoding guide-blocking and SapI changes are deliberately withheld for manual review, as are coding cases for which no verified synonymous solution is found.

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

## Important limitations

This is a computational design aid, not an experimentally validated clinical or diagnostic tool. Independently verify transcript choice, all junctions, the current plasmid sequence, sample genotype, guide activity, and biological suitability before ordering or experimentation. C-terminal tagging of Tubb5 may perturb its functionally important tubulin tail; the bundled result is computational only.

## Method and sequence sources

- Bollen et al. (2022), paper and supplementary S1-S3: guide priorities, 14-nt retargeting cutoff, 600-bp arms, SapI adapters, and mNeonGreen payload.
- Uploaded Addgene plasmid #169227 SnapGene file: fixed circular backbone sequence and annotation.
- Ensembl REST: live human/mouse transcript and reference-sequence retrieval.
- UCSC/GENCODE mm39: bundled Tubb5 reference fixture.
