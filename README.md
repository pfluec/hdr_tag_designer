# HDR Tag Designer

A local Streamlit prototype for Bollen-style **in-trans paired nicking (ITPN)** gene tagging with **SpCas9 D10A**. The first input is species: mouse **GRCm39** or human **GRCh38**.

Version 0.3.1 deliberately supports one finalized donor architecture:

- C-terminal `GGGGSAS-mNeonGreen-stop`
- Bollen TVBB C-term-mNeongreen
- Addgene plasmid **#169227**
- the supplied Addgene `.dna` sequence as the fixed backbone reference
- 600-bp arms by default
- no off-target analysis
- reference sequence only

Custom donor backbones and custom tags are intentionally deferred until after validation of this fixed-backbone test.

## What it does

1. Resolves a gene and transcript through public Ensembl REST in live mode.
2. Defines an N- or C-terminal insertion boundary from the selected transcript.
3. Finds all SpCas9-NGG candidates with a nominal nick within the chosen window.
4. Ranks guides by the Bollen priorities available locally: distance first, target disruption second, then basic GC/poly-T properties.
5. Generates gene-oriented homology arms and detects internal SapI sites.
6. For the bundled Tubb5 validation, applies one explicitly validated synonymous SapI-domestication change. No guide-blocking edit is added because the intended stop-codon replacement destroys the selected guide's PAM.
7. Generates the exact Bollen S1/S3 C-terminal UHA and DHA synthesis fragments and PCR-primer tail templates.
8. Parses the uploaded Addgene #169227 SnapGene file, verifies its four SapI sites and `TAC -> GGC -> TGA -> AAT` overhang order, and checks that its 729-bp payload matches Bollen supplementary S2.
9. Reconstructs the full circular Golden Gate product and verifies all four ligation junctions and the absence of residual SapI sites.
10. Exports TXT, CSV, FASTA, JSON, and an annotated GenBank file that can be opened in SnapGene and similar sequence editors.

The prototype does **not** calculate a validated on-target activity score, perform off-target searches, inspect sample-specific variants, or design locus-specific PCR annealing regions. Live designs that require a synonymous guide-blocking or SapI-domestication change are flagged for review rather than altered automatically; the automated corrections are currently validated only for the bundled Tubb5 test.

## Run locally

```zsh
cd hdr_tag_designer
./run.zsh
```

Or manually:

```zsh
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.txt
streamlit run app.py
```

## Run the tests and regenerate the bundled result

```zsh
./run_tests.zsh
```

Generated Tubb5 files are written to `outputs/tubb5_test/`.

## Bundled Tubb5 validation

The offline test uses mouse **Tubb5-201 / ENSMUST00000001566.10** on GRCm39/mm39. It validates the 444-aa coding endpoint and the terminal reference-genome interval used for the design.

The selected C-terminal guide is:

```text
5'-GAGGCAGAAGAGGAGGCCTA-AGG-3'
```

Its nominal D10A nick lies 1 bp from the insertion boundary. Stop-codon replacement deletes part of the guide's `AGG` PAM, so the complete target is absent from the edited locus and no additional guide-blocking mutation is required. One synonymous `GAG -> GAA` change removes an internal SapI site from the UHA.

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
