# HDR Tag Designer — Development Handoff

**Project version reviewed:** `0.6.0` (working tree; not yet committed)
**Handoff date:** 2026-07-21
**Purpose:** Continue development, debugging, and validation of the local HDR-tagging design tool in a new session or local coding environment.

---

## Future features — prioritized backlog

F1 was implemented in version `0.6.0`. F3 has a working baseline in `0.6.0`; its remaining genome-wide uniqueness check is documented below. F2 remains the next separate future feature.

### F1. Show the reference and final edited guide-binding sequence — IMPLEMENTED IN 0.6.0

Extend the **Selected guide** section so a user can directly compare what Cas9 sees before and after the donor edit.

Required presentation:

- reference target/protospacer and PAM, in guide orientation (`5' -> 3'`);
- final donor-allele sequence across the same genomic region after the intended insertion/deletion and every SapI/guide-blocking mutation;
- original and final PAM called out separately;
- strand and genomic coordinates;
- a compact alignment highlighting substitutions, deleted bases, the insertion boundary, and PAM changes;
- an explicit result such as `PAM destroyed`, `target split by insertion`, or the final longest retained target segment;
- the same information in TXT and JSON exports, not only in Streamlit.

An insertion can split the original 23-nt target and place hundreds of payload bases between its retained halves. Do not display a misleading synthetic 23-nt "final target" in that case. Show the actual edited-locus context, using a clearly labelled insertion marker/length in the compact alignment and an expandable full sequence where useful.

Implementation note: `_donor_target_after_edits()` and the retained-segment calculations already reconstruct part of this state internally. Promote the original/final target context into explicit result-model fields and test plus/minus strands, PAM destruction, synonymous PAM/seed changes, and targets split by N- or C-terminal insertion.

### F2. Twist Bioscience synthesis-complexity screening

Automatically screen every sequence intended for synthetic ordering, especially the final UHA and DHA synthesis fragments and any custom payload/orderable fragment, before labelling the design sequence-complete.

Required behavior:

- evaluate the final post-mutation sequences rather than the reference arms;
- identify the exact problematic interval and reason, rather than returning only a generic pass/fail;
- cover synthesis-relevant features such as extreme local/global GC, long homopolymers, direct or inverted repeats, low-complexity sequence, strong secondary-structure/hairpin potential, and other motifs/constraints in the current Twist Bioscience guidance;
- summarize each orderable fragment as `PASS`, `WARNING`, or `BLOCKED` and show the findings in Streamlit, TXT, and JSON;
- withhold an order-ready/sequence-complete label when a hard synthesis constraint is violated;
- where safe, propose alternative synonymous or non-coding mutations, but never silently change a coding sequence, splice-protected base, junction overhang, guide-protection mutation, or payload without rerunning every existing validation gate;
- retain the vendor/rule-set version and screening date in provenance.

Do not guess or permanently hard-code Twist thresholds from memory. At implementation time, verify the current official Twist synthesis guidelines and any available screening/API workflow, record the source/version, and keep vendor-specific thresholds configurable. Add deterministic fixtures for every rule so normal vendor/network availability is not required by the test suite.

Possible later extension: expose multiple eligible sequence-complexity/SapI fixes for user selection. Until an interactive mutation-choice workflow exists, keep safe changes automatic and make the exact selection rationale visible.

### F3. Genotyping primer design — BASELINE IMPLEMENTED IN 0.6.0

Design and report three genotyping PCR assays for every sequence-complete donor design:

1. **Non-edited/WT locus assay**
   - Design primers that amplify the non-edited locus across or around the native insertion site.
   - Prefer a configuration that clearly distinguishes the WT allele from the edited allele by allele specificity or expected amplicon size.
   - Report whether the edited allele is expected to amplify, and give the expected WT and edited product sizes where applicable.
2. **5-prime insertion-junction assay**
   - Use a locus-specific forward primer in genomic sequence on the 5-prime side, preferably outside the UHA so that residual donor plasmid or correctly assembled donor DNA alone cannot produce the expected product.
   - Use a reverse primer inside the inserted payload.
   - Because the reverse primer depends only on the payload/backbone architecture, define and validate it once per payload configuration and reuse it across genes using the same targeting backbone where the junction geometry permits.
3. **3-prime insertion-junction assay**
   - Use a forward primer inside the inserted payload.
   - Use a locus-specific reverse primer in genomic sequence on the 3-prime side, preferably outside the DHA.
   - Because the forward primer depends only on the payload/backbone architecture, define and validate it once per payload configuration and reuse it across genes using the same targeting backbone where the junction geometry permits.

Required output and validation:

- primer names, sequences (`5' -> 3'`), strand/orientation, genomic or payload coordinates, length, GC percentage, melting temperature, and expected amplicon sequence/size;
- explicit indication of which primers are locus-specific and which are reusable payload primers;
- separate predictions for WT and edited alleles, including the reason an assay should or should not amplify each allele;
- specificity checks against the selected reference genome and against the complete assembled donor plasmid;
- rejection or warning for problematic homopolymers, repeats, strong self-dimers, cross-dimers, hairpins, extreme GC/Tm, or non-unique genomic binding;
- primer placement outside the homology arms for external junction primers whenever sequence quality and amplicon-size constraints allow;
- plus- and minus-strand genes handled without changing the biological 5-prime/3-prime assay definitions;
- Streamlit display plus TXT, CSV, and JSON exports;
- configurable target amplicon sizes and PCR/primer constraints, with the rule-set and genome assembly recorded in provenance.

Implementation should use the final edited locus and assembled payload sequences, including every automatic SapI and guide-blocking mutation. Reusable payload primers must be keyed to the actual payload definition/version rather than only to a human-readable backbone name. Add regression fixtures for N- and C-terminal designs, both gene strands, WT discrimination, repeated/non-unique genomic primer sites, and custom payloads.

---

## 0. Continuation status — read this first

Development is continuing on the Git branch:

```text
feat/custom-payloads-guide-context-genotyping
```

Version `0.6.0` retains the automatic-mutation and N/C backbone work from `0.5.1` and adds:

1. complex custom payload cassettes without a forced single-ORF interpretation;
2. explicit warnings for non-frame, multi-ORF, missing-linker, or internal-stop payloads;
3. reference, post-point-mutation, compact actual-edited, and full actual-edited guide-target sequences;
4. WT-locus, 5-prime junction, and 3-prime junction primer design with Primer3 thermodynamics;
5. junction genomic primers outside the homology arms and payload primers at least 150 bp from the tested junction;
6. genotyping display plus TXT, CSV, FASTA amplicon, and JSON exports.
7. annotated WT and edited linear locus contexts extending 300 bp beyond both homology arms, with arm, payload/native-junction, and genotyping-primer features.
8. download buttons use Streamlit's no-rerun mode, and the latest completed `DesignResult` is restored from session state after any incidental rerun instead of resetting the page.

The previous release's behavior remains intact:

- bounded Ensembl retries for HTTP 429, HTTP 5xx, and connection failures;
- clean user-facing Ensembl errors without returned HTML bodies;
- validation that an explicit transcript belongs to the selected gene and is protein-coding;
- a user decision encoded as a regression rule: **guide proximity remains the primary ranking criterion, and a nearer guide is not demoted merely because it needs a verified silent mutation**;
- a SapI quality-control panel with total/per-arm counts and site-level mutation details;
- mutation propagation into final arms, edited CDS, plasmid assembly, TXT, JSON, FASTA, and GenBank output;
- a conda-aware test runner and `pytest` in the environment/test dependency specifications.

Current validation state:

```text
36/36 unittest tests pass
Bundled Tubb5 fixture remains sequence-complete
Exact #169227 Tubb5 sequence/coordinate/plasmid regressions pass unchanged
Synthetic N-terminal designs pass on plus and minus gene strands
Addgene #169226 yields a 3947-bp SapI-free assembled plasmid with 600-bp arms
Custom-upload classification and complete assembly path pass using #169226 as a fixture
Complex multi-ORF/non-frame payload classification and non-blocking assembly are covered
Tubb5 produces all three constrained genotyping primer pairs
Live Actb N-terminal validation remains sequence-complete and produces all three primer pairs
```

The bundled offline Tubb5 fixture remains `ENSMUST00000001566.10` for reproducibility. Ensembl resolved the same stable mouse transcript as version `.11` during the 2026-07-21 live check.

### 0.1 Automatic mutation behavior now implemented

For the selected nearest guide, the live design path:

1. re-evaluates any existing SapI-domestication changes against the guide;
2. prefers synonymous destruction of either invariant `G` in the NGG PAM;
3. if PAM destruction is impossible, searches synonymous protospacer/seed changes that reduce the longest retained segment to `<=14 nt`;
4. searches combinations affecting up to three codons;
5. ranks non-PAM solutions by fewest changed bases, proximity to the PAM, then remaining segment length;
6. reconstructs and translates the complete CDS and requires an identical protein sequence;
7. rejects changes that introduce a SapI site, lengthen a homopolymer to at least six bases, or touch the first/last three exonic bases near a splice boundary;
8. leaves noncoding or otherwise unresolved **guide-blocking** changes blocked instead of silently switching to a farther guide.

Generic SapI domestication examines every `GCTCTTC` and `GAAGAGC` motif in both homology arms. It first applies the smallest verified synonymous coding change available. If none exists, it searches one-base non-coding substitutions, preferring bases outside the mature transcript over UTR bases and excluding every coding base, the first/last three bases of each exon, and six intronic bases on either side of each splice boundary. A released candidate must remove the original site without creating a new SapI site or lengthening a problematic homopolymer. The fixed Tubb5 fixture retains its exact locked `GAG -> GAA` regression correction.

Live regression check for mouse **Actb**, N-terminal mode, canonical `ENSMUST00000100497.11`:

- original UHA SapI site: `GCTCTTC`, arm bases 155-161, `chr5:142,891,878-142,891,884`;
- location: intron between the first two transcript exons;
- automatic change: UHA base 158, `chr5:142,891,881`, gene-oriented `C>A`;
- final arm SapI sites: zero;
- final result: sequence-complete 3,947-bp N-terminal plasmid with zero residual SapI sites.

The former generic `REVIEW REQUIRED` status was renamed to `DESIGN BLOCKED - NO SEQUENCE-COMPLETE OUTPUT`, because the current UI has no sequence-editing/review action. A future mutation-choice interface may expose the eligible non-coding SapI candidates; version 0.6.0 selects automatically and records the exact ranking rationale.

### 0.2 N-terminal and custom-backbone behavior now implemented

The newly supplied `data/bollen_supplementary_s1.docx` and `data/addgene-169226.dna` establish the N-terminal architecture:

- insert immediately after the endogenous `ATG`;
- UHA adapters produce `TAC -> GTG`;
- DHA adapters produce `AGC -> AAT`;
- #169226 is 2,765 bp, circular, SHA-256 `75d9d25b4dac8083c401ee5ac76b080a5f62e42d23357a81b4ac33e84a434177`;
- the supplied protocol file SHA-256 is `5270e6028876d74daf551e4607d274ce0411fc07e1caaead3071bdb9378f44e0`;
- inner cuts reconstruct a 726-nt payload: 705-nt mNeonGreen followed by the 21-nt GGGGSAS linker, with no stop codon;
- a standard 600/600-bp design produces a 1,972-bp donor insert and 3,947-bp final plasmid;
- predicted edited protein order is native initiator methionine, payload, then native residue 2 onward.

Custom `.dna` uploads currently require all of the following structural properties before release:

- circular SnapGene input;
- exactly four SapI sites;
- exact supported order `TAC/GGC/TGA/AAT` (C) or `TAC/GTG/AGC/AAT` (N);
- cassette coordinates in increasing file order (a cassette crossing the sequence origin must be rotated first);
- selected UI terminus matching the inferred backbone terminus.

Payload frame, linker, ambiguity, and internal-stop checks are now interpretation checks rather than release gates. A conventional in-frame cassette with the architecture-specific GGGGSAS linker receives a fusion translation. Every other cassette is labelled `complex cassette`, produces a visible warning, and is assembled without asserting a single fusion protein. Passing uploads use the same full donor/plasmid assembly, junction checks, residual-SapI check, exports, and QC gates as the built-ins.

### 0.3 Guide-context and genotyping behavior added in 0.6.0

Every ranked guide now records the reference 23-nt target, the same target after all donor point mutations, and the actual donor-allele sequence across the original genomic span. When an insertion splits the target, the compact representation includes `[INSERT n nt]`; the full inserted sequence is retained separately. This prevents presenting a misleading synthetic 23-nt final target. Tubb5 C-terminal and synthetic N-terminal plus/minus-strand regressions cover the reconstruction. The live Actb result shows its synonymous `GAC -> GAT` guide-blocking change directly in the post-mutation target.

Genotyping primer design is implemented in `hdr_designer/primers.py` with Primer3/primer3-py. The current fixed rules are 18-27 nt, 57-63 C, 35-65% GC, no homopolymer longer than four, pair Tm difference at most 3 C, and hairpin/homodimer/heterodimer Tm below 45 C. Candidate pairs must be unique within the constructed assay template. External locus primers bind in a 300-bp fetched region outside the relevant arm and must be absent in both orientations from the complete assembled donor plasmid. Payload primers are keyed by SHA-256 of the actual payload and prioritized deterministically for reuse across loci where a compatible genomic partner exists.

The WT assay uses two external genomic primers and distinguishes alleles by product size. The 5-prime junction assay uses an external 5-prime genomic forward primer plus a payload reverse primer whose binding start is at least 150 bp from the 5-prime junction. The 3-prime assay uses a payload forward primer whose binding end is at least 150 bp from the 3-prime junction plus an external 3-prime genomic reverse primer.

Important remaining limitation: genome-wide primer uniqueness is not computed locally. The app marks this check `NOT RUN` and tells users to confirm every pair with Primer-BLAST before ordering. Donor-plasmid specificity and uniqueness within the modeled assay template are enforced automatically. This limitation should not be silently removed from the UI or reports until a real reference-genome specificity backend exists.

### 0.4 Important implementation locations

- `hdr_designer/backbones.py`: built-in definitions, SapI classification, payload extraction, custom `.dna` inference, shared synthesis fragments, and shared plasmid assembly.
- `hdr_designer/design.py`: transcript/genome mapping, arm generation, generic synonymous mutation search, guide blocking, SapI domestication, N/C fusion placement, and final validation gates.
- `hdr_designer/primers.py`: constrained WT/junction primer search, Primer3 thermodynamics, assay geometry, payload reuse keys, and donor-plasmid checks.
- `hdr_designer/ensembl.py`: bounded retry behavior and transcript-parent/biotype validation.
- `hdr_designer/exports.py`: mutation-aware exports, guide context, genotyping-primer CSV/report/FASTA output, and `sapi_qc_rows()`.
- `app.py`: dynamic N/C built-in selection, custom `.dna` upload, SapI QC, guide-context comparison, and genotyping-primer display.
- `tests/test_online_design.py`: PAM blocking, seed blocking, unsafe refusal, one/multiple SapI sites, and reverse-strand coverage.
- `tests/test_guides.py`: nearest-guide precedence regression.
- `tests/test_ensembl.py`: retry and transcript-validation coverage.
- `tests/test_tubb5.py`: fixed sequence/coordinate, assembly, export, SapI-QC, guide-context, and genotyping-primer regression.

### 0.5 Next development objective

The next chunk of work is now:

> Separate payload metadata from backbone configuration, then add payload-only DNA/FASTA input and stronger negative fixtures without regressing either built-in architecture.

Recommended sequence:

1. introduce a dedicated serializable `PayloadDefinition` rather than storing tag/linker offsets directly on `BackboneDefinition`;
2. support payload-only DNA/FASTA input for either N- or C-terminal architecture;
3. allow GenBank input in addition to SnapGene `.dna`;
4. add negative fixtures for malformed files, wrong site count/order, residual SapI, frame error, internal stop, missing linker, terminus mismatch, and origin-crossing cassettes;
5. decide whether origin-crossing circular cassettes should be normalized automatically;
6. add user-provided payload/tag names and export annotations;
7. keep the exact #169227 Tubb5 and #169226 N-terminal regression outputs locked.

Sections later in this document describe the original `0.4.0` roadmap and should be read as historical design context where they conflict with this continuation-status section.

---

## 1. Project objective

The goal is to build a local Streamlit application that designs CRISPR–Cas9 HDR donor constructs for N- or C-terminal endogenous gene tagging in mouse and human cells.

The initial implementation follows the workflow described by **Bollen et al. (2022)**, especially the supplementary protocol describing:

- guide selection near the intended insertion site,
- homology-arm design,
- placement of guide target sequences outside the homology arms for in-trans paired nicking,
- SapI-compatible Golden Gate assembly,
- sequence validation of the donor and edited locus.

The current prototype has two validated built-in donor architectures and a deliberately narrow custom-backbone path before broader user-defined payload support.

---

## 2. Confirmed user requirements

The following choices were explicitly agreed upon:

- **Species:** mouse and human
- **Reference assemblies:** mouse GRCm39 and human GRCh38
- **Editing strategy:** in-trans paired nicking using SpCas9 D10A
- **Default arm length:** 600 bp
- **Guide nuclease motif:** SpCas9 NGG
- **Off-target analysis:** not required at this stage
- **Sequence source:** public Ensembl reference sequence
- **Sample-specific variants:** not included
- **Interface:** local Streamlit application
- **Built-in donors:** Addgene plasmids **#169226** (N-terminal) and **#169227** (C-terminal)
- **Built-in tag:** mNeonGreen with architecture-specific GGGGSAS linker placement
- **Validation gene:** mouse **Tubb5**
- **Validation transcript:** `ENSMUST00000001566`
- **Custom backbones:** supported for structurally matching circular SnapGene `.dna` files
- **Custom payload-only input:** next development chunk
- **Blocking mutations:** automatically design verified synonymous coding changes when possible; otherwise block sequence-complete output

---

## 3. Current package

The current repository version is:

```text
0.6.0 working tree on feat/custom-payloads-guide-context-genotyping
```

Main project structure:

```text
hdr_tag_designer/
├── app.py
├── requirements.txt
├── pyproject.toml
├── run.zsh
├── run_tests.zsh
├── README.md
├── VALIDATION.txt
├── VALIDATION_SUMMARY.txt
├── hdr_designer/
│   ├── __init__.py
│   ├── models.py
│   ├── sequence.py
│   ├── ensembl.py
│   ├── guides.py
│   ├── design.py
│   ├── backbones.py
│   ├── snapgene.py
│   ├── exports.py
│   └── fixtures.py
├── data/
│   ├── addgene_169227.dna
│   ├── bollen_169227_c_terminal_payload.fa
│   ├── bollen_s2_mneongreen_c_terminal_fragment.txt
│   └── tubb5_ucsc_mrna.fa
├── scripts/
│   └── generate_tubb5_report.py
└── tests/
    ├── test_tubb5.py
    ├── test_online_design.py
    ├── test_ensembl.py
    ├── test_guides.py
    └── test_app.py
```

---

## 4. Current functional scope

### 4.1 Supported inputs

The Streamlit app currently accepts:

1. species,
2. gene symbol or Ensembl gene ID,
3. optional explicit Ensembl transcript ID,
4. N- or C-terminal mode,
5. fixed donor display,
6. homology-arm length,
7. guide-search radius.

For reliable use, explicitly supplying the desired transcript ID is currently recommended.

### 4.2 Sequence retrieval

The Ensembl client currently:

- resolves gene symbols through `/xrefs/symbol/`,
- retrieves expanded gene and transcript records through `/lookup/id/`,
- retrieves cDNA and CDS sequences through `/sequence/id/`,
- retrieves genomic regions through `/sequence/region/`,
- strips transcript version suffixes before lookup,
- supports mouse and human.

The explicit transcript workflow appears to work more reliably than allowing the application to resolve the transcript automatically.

### 4.3 Terminus definition

For C-terminal tagging:

- the insertion is placed immediately before the endogenous stop codon,
- the endogenous stop codon is removed from the downstream arm,
- the fixed donor payload supplies its own stop codon.

For N-terminal tagging:

- the insertion boundary is placed immediately after the start codon,
- the current version only provides a locus/guide/arm preview,
- complete plasmid assembly is not validated because Addgene #169227 is a C-terminal backbone.

### 4.4 Guide discovery

The current guide engine:

- enumerates SpCas9 NGG targets on both strands,
- calculates a nominal D10A nick boundary,
- searches within a user-defined distance of the insertion site,
- ranks candidates by:
  1. nick distance from insertion,
  2. whether the intended edit disrupts the guide,
  3. poly-T avoidance,
  4. GC content close to 50%.

The current version does **not** run:

- a validated on-target activity model,
- genome-wide off-target analysis,
- bulge-aware off-target analysis.

The guide scoring in the interface should therefore be treated as a transparent heuristic, not a full guide-quality prediction.

---

## 5. Bollen guide-retargeting rule currently represented

The code uses a guide safety cutoff of:

```text
14 nt
```

The intended logic is:

- if the intended edit destroys the PAM, no additional blocking mutation is required;
- if the PAM remains intact, calculate the longest uninterrupted retained segment of the original target;
- if more than 14 nt remain contiguous with an intact PAM, the donor must contain an additional guide-blocking change;
- if the design requires a blocking mutation, attempt a verified synonymous PAM or seed change;
- if no safe synonymous coding change is available, withhold final order-ready cloning fragments.

Version `0.4.0` implements this for coding portions of live mouse and human designs. Noncoding changes remain manual-review cases. The selected nearest guide is retained while a silent blocking solution is sought; a farther guide is not silently substituted.

---

## 6. Homology-arm logic

The current version generates gene-oriented 5′ and 3′ homology arms.

Default:

```text
600 bp UHA
600 bp DHA
```

For reverse-strand genes, chromosome-forward genomic sequence is reverse-complemented so all displayed arm and donor sequences are in gene orientation.

The tool stores both:

- chromosome-forward reference sequence,
- gene-oriented sequence.

The interface should continue to make this distinction explicit because reverse-strand loci are a major source of design errors.

---

## 7. Fixed Bollen donor architecture

The current fully assembled design is locked to:

```text
Addgene #169227
C-terminal GGGGSAS-mNeonGreen-stop
```

The supplied SnapGene `.dna` file is bundled as:

```text
data/addgene_169227.dna
```

Previously verified fixed-backbone properties:

```text
Backbone length: 2,768 bp
Topology: circular
SapI sites: 4
Expected overhang order: TAC -> GGC -> TGA -> AAT
Payload length: 729 bp
Final payload: GGGGSAS linker + 235-aa mNeonGreen + TGA
```

The SHA-256 checksum recorded for the supplied backbone was:

```text
b5ccaa5a257b71a1f2bac05ab15785f07098a46ed805c2862b5beeade04046b1
```

---

## 8. SapI synthesis-fragment architecture

The SapI extensions are taken directly from the Bollen supplementary protocol.

For the current C-terminal donor system:

### UHA synthesis fragment

```text
AACGCTCTTCATAC
+ target-with-PAM
+ UHA
+ GGCTGAAGAGCGCG
```

### DHA synthesis fragment

```text
CGCGCTCTTCGTGA
+ DHA
+ target-with-PAM
+ AATCGAAGAGCGTT
```

Expected SapI-generated overhangs:

```text
UHA fragment: TAC / GGC
DHA fragment: TGA / AAT
```

The current program verifies the overhangs and simulates the final Golden Gate product.

---

## 9. Tubb5 computational validation

The bundled validation case is:

```text
Species: mouse
Gene: Tubb5
Transcript stable ID: ENSMUST00000001566
Fixture version: .10
Assembly: GRCm39/mm39
Terminus: C-terminal
Arm length: 600 bp
```

### Selected guide

```text
Spacer: GAGGCAGAAGAGGAGGCCTA
PAM: AGG
Target + PAM: GAGGCAGAAGAGGAGGCCTAAGG
PAM-containing chromosome strand: -
Nick distance from insertion: 1 bp
```

The intended replacement of the endogenous stop codon destroys part of the PAM, so no additional guide-blocking mutation is required for this test case.

### Tubb5 arms

```text
UHA: chr17:36,145,877–36,146,476
DHA: chr17:36,145,274–36,145,873
```

Coordinates above are represented as 1-based intervals in the validation report.

### SapI domestication

The Tubb5 UHA contains an internal SapI site.

The bundled fixture applies:

```text
GAG -> GAA
Glu -> Glu
```

This is a synonymous change and removes the internal SapI site.

The exact fixture correction remains locked as a regression test. Live designs now use the generic synonymous SapI domesticator for coding arm sites; noncoding or unresolved sites remain blocked.

### Predicted fixed fusion

```text
Tubb5: 444 aa
Linker: 7 aa, GGGGSAS
mNeonGreen: 235 aa
Predicted fusion: 686 aa
```

### Simulated assembled plasmid

```text
Donor insert: 1,975 bp
Final circular plasmid: 3,950 bp
Residual SapI sites: 0
```

This result is a computational test only. C-terminal tagging may interfere with the functional tubulin tail.

---

## 10. Current exports

The program currently exports combinations of:

- plain-text design report,
- guide table as CSV,
- homology arms and cloning fragments as FASTA,
- JSON design data,
- annotated GenBank assembled plasmid,
- junction sequences,
- complete homology-arm cloning primers with separately reported SapI/Golden Gate tails and genomic annealing regions.

The Streamlit interface and TXT report also include SapI quality control: total/per-arm counts, every original motif and coordinate, resolution status, exact nucleotide/codon change, protein consequence, and selection reason.

Version 0.6.0 appends Primer3-evaluated, arm-end genomic annealing regions to all four fixed Bollen cloning tails. The UI colors tails and annealing regions separately and provides a copyable complete sequence. If a required donor mutation lies internally, outside both endpoint-primer footprints, the app warns that endpoint PCR alone cannot create the final arm and that synthesis or an additional mutagenesis strategy is required.

---

## 11. Environment and installation

The project now uses the named conda environment `hdr-tag-designer`. Both `run.zsh` and `run_tests.zsh` activate that environment rather than installing into an unrelated active Python environment. `environment.yml` specifies Python 3.11 and includes the runtime dependencies plus pytest.

### Recommended local environment

Create an `environment.yml`:

```yaml
name: hdr-tag-designer

channels:
  - conda-forge

dependencies:
  - python=3.11
  - pip
  - pip:
      - streamlit>=1.59,<2
      - pandas>=2.0,<3
      - requests>=2.31,<3
      - biopython>=1.83,<2
      - pytest
```

Create and activate it:

```zsh
conda env create -f environment.yml
conda activate hdr-tag-designer
streamlit run app.py
```

For iterative development:

```zsh
python -m pip install -e .
pytest
streamlit run app.py
```

The already-created local environment did not contain pytest during this continuation, so the authoritative suite was run through the repository's `unittest` runner. `./run_tests.zsh` passes all 28 tests and regenerates the Tubb5 outputs.

---

## 12. Ensembl transcript resolution

An observed failure was:

```text
Ensembl REST returned HTTP 500 for /lookup/id/ENSMUST...
```

The app worked when the desired transcript ID was supplied explicitly.

Version `0.4.0` adds bounded retries for HTTP 429, HTTP 5xx, and connection failures, plus concise errors that do not expose returned HTML. Explicit transcript IDs are checked against the selected parent gene and protein-coding biotype.

Still to implement:

1. display all current protein-coding transcripts for user selection;
2. use the canonical transcript only as a visible default, not as a hidden assumption;
3. show transcript version, protein length, CDS length, and terminal exon;
4. add an archive-ID fallback for retired Ensembl IDs;
5. cache successful Ensembl responses locally;
6. allow a locally supplied sequence/GenBank record when Ensembl is unavailable.

---

# 13. Implemented feature: automatic guide-blocking mutations

## 13.1 Goal

When the intended tag insertion does not sufficiently destroy the selected guide target, the program should alter the donor homology arm so the correctly edited allele is resistant to further nicking.

The mutation must not alter the desired protein or disrupt essential local sequence features.

**Status in 0.4.0:** implemented for synonymous coding changes in the selected guide. PAM destruction is preferred; otherwise the engine searches PAM-proximal protospacer changes that satisfy the 14-nt cutoff. Complete-CDS translation, SapI avoidance, a three-base exonic splice-edge exclusion, and homopolymer checks are enforced. Codon-usage scoring and comprehensive cryptic-splice/regulatory-motif prediction are not yet implemented. Noncoding changes are never released automatically.

## 13.2 Suggested mutation priority

For each candidate guide requiring blocking:

1. **Destroy the PAM**
   - Prefer a synonymous coding change if the PAM overlaps CDS.
   - In noncoding sequence, change the PAM directly where biologically acceptable.

2. **Mutate PAM-proximal seed bases**
   - Prefer positions closest to the PAM.
   - Use the minimum number of mutations required.
   - Prefer synonymous changes in coding sequence.

3. **Break the longest retained target segment**
   - Recalculate the longest contiguous original target sequence after every proposed change.
   - Require the final retained segment to satisfy the Bollen cutoff.

4. **Re-evaluate the complete edited locus**
   - Confirm that a functional NGG target is no longer present.
   - Confirm that the donor target sites outside the arms remain as intended.

## 13.3 Coding-sequence mutation algorithm

For mutations within CDS:

1. map each genomic guide/PAM base to transcript and CDS position;
2. identify the overlapping codon;
3. enumerate synonymous codons for that amino acid;
4. score possible codon changes by:
   - PAM destruction,
   - seed disruption,
   - number of nucleotide changes,
   - codon usage,
   - avoidance of new SapI sites,
   - avoidance of long homopolymers,
   - avoidance of cryptic splice motifs near exon boundaries;
5. reconstruct the complete CDS;
6. translate both original and edited CDS;
7. require identical protein sequence apart from the intended tag.

## 13.4 Noncoding mutation handling

If the guide overlaps:

- intron,
- UTR,
- intergenic sequence,
- untranslated part of the terminal exon,

the program should not automatically assume that any nucleotide can be changed safely.

Add configurable exclusion zones for:

- splice donor and acceptor motifs,
- branch-point/polypyrimidine regions where annotated,
- Kozak sequence near the start codon,
- stop-codon context,
- known regulatory motifs when annotations are available.

At minimum, noncoding blocking changes should be marked for manual review.

## 13.5 Required output

The final report should explicitly show:

```text
Mutation purpose: guide blocking
Genomic coordinate
Homology arm and arm-relative position
Original base/codon
Edited base/codon
Protein consequence
PAM before and after
Longest retained target segment before and after
Reason this mutation was selected
```

The arm FASTA, donor plasmid, edited-locus simulation, and GenBank annotations must all include the blocking mutation.

---

# 14. Implemented feature: generic SapI domestication

Version `0.4.0` automatically fixes verified synonymous coding SapI sites in generic live designs. The remaining work in this section applies to custom payloads/backbones, noncoding sites, and broader type-IIS enzyme support.

The general implementation should:

1. identify every `GCTCTTC` and `GAAGAGC` site in:
   - UHA,
   - DHA,
   - custom payload,
   - final assembled donor;

2. map the site to genomic/transcript/CDS coordinates;

3. when coding:
   - enumerate synonymous substitutions;
   - choose the smallest safe edit;
   - preserve protein sequence;

4. when noncoding:
   - propose candidate edits but require review;

5. check that each correction:
   - removes the original SapI site,
   - does not create another SapI site,
   - does not recreate a functional guide target,
   - does not create an unwanted type-IIS site if broader enzyme screening is added;

6. rerun full plasmid assembly and translation validation.

If a safe mutation cannot be found, possible fallbacks are:

- shift the homology-arm boundary,
- shorten or lengthen the arm,
- use a different guide,
- use a different type-IIS assembly strategy,
- use a synthesized donor with a different cloning approach.

---

# 15. Priority feature: custom donor plasmids

## 15.0 Current starting point and recommended phases

The fixed implementation in `hdr_designer/backbones.py` already provides a validated reference behavior:

- SnapGene parsing through `hdr_designer/snapgene.py`;
- exact Addgene #169227 checksum, topology, length, four SapI sites, and overhang order;
- payload extraction and comparison to the Bollen S2 fixture;
- full circular assembly and feature-coordinate reconstruction;
- final SapI and junction validation.

Do not replace this in one step. Convert it into the first `BackboneDefinition` configuration while preserving the exact sequences, coordinates, 3,950-bp final plasmid, and GenBank annotations asserted by `tests/test_tubb5.py`.

Recommended phases:

1. **Configuration models without UI changes**
   - Add `BackboneDefinition`, `PayloadDefinition`, `AssemblyJunction`, and validation-result models.
   - Use dataclasses and JSON-serializable dictionaries first; YAML can be added later without making PyYAML a prerequisite.
   - Express Addgene #169227 as the built-in C-terminal configuration and route the existing fixed workflow through it.

2. **Normalized uploaded sequence record**
   - Define one internal sequence/feature/topology representation.
   - Adapt the existing SnapGene parser to uploaded bytes or a safely managed temporary file.
   - Parse GenBank with Biopython and FASTA with intentionally reduced annotation guarantees.
   - Enforce file-size limits, DNA validation, topology choice, and actionable parse errors.

3. **User-specified assembly configuration before automatic inference**
   - Let the user identify the enzyme, cassette boundaries/sites, expected overhangs, payload interval/orientation, terminus, and target-site placement.
   - Validate those declarations against the uploaded sequence.
   - Only after this path is reliable should the app attempt automatic cassette inference.

4. **General assembly simulator**
   - Parameterize recognition motif, cut offsets, overhang length/order, retained vector intervals, and feature shifts.
   - Start with SapI as the only released enzyme configuration; supporting arbitrary type-IIS enzymes requires a validated enzyme-definition model and dedicated fixtures.

5. **Release gates and UI**
   - Never produce order-ready fragments unless topology, site count, overhang uniqueness/order, payload orientation/frame, internal-site checks, circular reconstruction, and junction checks all pass.
   - Show a backbone QC panel analogous to the new SapI arm QC panel.

Suggested first acceptance test:

> Loading the bundled Addgene #169227 file through the new generic configuration path must produce byte-for-byte identical synthesis fragments, donor insert, assembled plasmid sequence, and equivalent GenBank annotations to the current fixed path.

## 15.1 Desired workflow

Allow the user to upload a donor backbone in one of these formats:

- SnapGene `.dna`
- GenBank `.gb` or `.gbk`
- FASTA, with reduced annotation support

The user should then identify:

- donor topology,
- type-IIS enzyme,
- left and right assembly sites,
- payload region,
- linker,
- tag CDS,
- stop-codon behavior,
- selection cassette if present,
- N- or C-terminal architecture,
- expected assembly overhangs,
- target-site placement outside the homology arms.

## 15.2 Automatic backbone analysis

The application should attempt to:

1. parse the sequence and annotations;
2. detect all SapI or selected type-IIS sites;
3. calculate the generated overhangs;
4. identify the replaceable cloning cassette;
5. verify that each overhang is unique and compatible;
6. reconstruct the unmodified payload;
7. classify payload reading-frame/linker compatibility without rejecting complex cassettes;
8. confirm circular assembly;
9. report unexpected internal restriction sites.

## 15.3 Backbone configuration object

Backbone logic should be moved from hard-coded constants into a serializable configuration such as:

```yaml
name: Bollen TVBB C-term mNeonGreen
addgene_id: 169227
topology: circular
terminus: C
enzyme: SapI
linker_sequence: GGGGSAS
payload_supplies_stop: true
overhangs:
  vector_to_uha: TAC
  uha_to_payload: GGC
  payload_to_dha: TGA
  dha_to_vector: AAT
target_site_position:
  uha: before_arm
  dha: after_arm
```

The exact schema can be refined during implementation.

## 15.4 Safety checks for custom backbones

Before producing order-ready sequences, require:

- expected number of type-IIS sites,
- unique assembly overhangs,
- no internal enzyme sites in final fragments,
- valid donor topology,
- correct payload orientation,
- explicit reading-frame/linker/stop interpretation,
- complete circular-plasmid reconstruction where applicable.

---

# 16. Priority feature: custom tags and payloads

## 16.0 Recommended first release scope

Implement custom payloads in a deliberately narrow order:

1. accept a named DNA cassette by paste or FASTA upload;
2. accept optional linker DNA and explicit single-fusion versus complex-cassette intent;
3. validate DNA alphabet, report frame/translation/internal stops, and enforce terminus, junction, and SapI compatibility without rejecting multi-ORF or non-frame cassettes;
4. pass a classified `PayloadDefinition` into the generalized backbone assembler;
5. propagate the payload name, checksum, sequence, translation, and validation decisions into JSON/TXT/FASTA/GenBank output.

Defer protein-to-DNA reverse translation and codon optimization until the user can choose a target organism/codon table and the project has exact reproducibility tests. GenBank-feature selection and a local payload library can follow the initial DNA/FASTA path. Selection cassettes, artificial introns, and P2A architectures should remain explicit advanced configurations rather than being inferred from an arbitrary sequence.

Allow the user to supply:

- DNA sequence,
- protein sequence with reverse translation,
- GenBank feature,
- FASTA sequence,
- tag selected from a local library.

Configuration should include:

- tag name,
- coding sequence,
- linker sequence,
- whether the tag includes start or stop codons,
- whether the tag is N- or C-terminal,
- optional flexible linker,
- optional cleavage peptide,
- optional selection cassette,
- optional artificial intron,
- codon-optimization setting,
- fluorescence/protein metadata.

Required checks:

- tag sequence length divisible by three,
- translation contains no unintended internal stop,
- correct frame at both junctions,
- no incompatible SapI sites,
- no duplicated start or stop codon,
- expected fusion-protein length,
- exact junction translation.

Useful preset payloads later could include:

- mNeonGreen,
- mScarlet,
- EGFP,
- mStayGold,
- HaloTag,
- SNAP-tag,
- FLAG,
- HA,
- V5,
- P2A-selection cassettes.

---

# 17. N-terminal full donor support

Current N-terminal mode is a preview only.

To make it complete:

1. obtain and validate the appropriate Bollen N-terminal TVBB backbone;
2. load the corresponding `.dna` or GenBank sequence;
3. implement the N-terminal SapI adapter arrangement;
4. determine whether the endogenous ATG is retained or supplied by the donor;
5. validate Kozak/start-codon context;
6. reconstruct the complete N-terminal fusion CDS;
7. handle transcripts with alternative start codons;
8. verify that the tag does not introduce an upstream start or frame shift;
9. add N-terminal fixture tests.

N-terminal and C-terminal architectures should be separate backbone configurations rather than conditional hard-coded strings.

---

# 18. Transcript selection improvements

A future transcript-selection panel should show:

- stable transcript ID and version,
- MANE status for human where applicable,
- Ensembl canonical status,
- protein-coding status,
- CDS length,
- protein length,
- terminal exon coordinates,
- start and stop genomic positions,
- whether alternative transcripts share the same insertion boundary,
- whether the selected tag would affect all or only some isoforms.

The user should be warned when:

- alternative transcripts use different final exons,
- the selected terminus is not shared,
- the transcript is noncoding,
- the CDS lacks a standard ATG or stop codon,
- the transcript ID is retired or no longer current.

---

# 19. Additional useful features

## 19.1 PCR primer design — implemented baseline

Add locus-specific annealing sequences for:

- UHA amplification,
- DHA amplification,
- junction PCR,
- external validation PCR,
- optional sequencing primers.

Use a primer-design library or Primer3 and report:

- melting temperature,
- GC content,
- amplicon size,
- product specificity where possible,
- primer-dimer/hairpin warnings.

## 19.2 Sequence and variant inputs

Support:

- pasted genomic sequence,
- uploaded FASTA,
- uploaded GenBank,
- uploaded VCF,
- selected mouse strain variants,
- cell-line sequence overrides.

Flag variants within:

- guide spacer,
- PAM,
- homology arms,
- primer annealing sites.

## 19.3 Guide scoring

Off-target analysis remains optional, but future guide ranking could include:

- Doench-style on-target score,
- CRISPRscan or another transparent model,
- poly-T and extreme-GC penalties,
- terminal nucleotide preferences,
- repetitive-sequence filters.

Keep guide proximity and target disruption as prominent Bollen-specific criteria.

## 19.4 Caching and reproducibility

Record with every design:

- species,
- assembly,
- Ensembl release or retrieval date,
- transcript stable ID and version,
- all genomic coordinates,
- backbone checksum,
- payload checksum,
- application version,
- configuration file,
- mutation decisions.

Cache Ensembl responses so a design can be reproduced if online annotations later change.

## 19.5 Improved error handling

Replace raw exceptions with actionable messages:

- no gene found,
- transcript does not belong to gene,
- temporary Ensembl failure,
- noncoding transcript,
- nonstandard start/stop codon,
- no guide in search window,
- unresolved SapI site,
- blocking mutation required,
- backbone structural mismatch.

---

# 20. Suggested refactoring

The prototype currently contains fixed assumptions spread through `design.py` and `backbones.py`.

A cleaner architecture would separate:

```text
SequenceProvider
TranscriptResolver
InsertionSiteStrategy
GuideFinder
GuideBlockingDesigner
HomologyArmDesigner
RestrictionSiteDomesticator
PayloadDefinition
BackboneDefinition
AssemblySimulator
DesignValidator
Exporter
```

Suggested data models:

```text
BackboneDefinition
PayloadDefinition
AssemblyJunction
SequenceMutation
GuideBlockingResult
RestrictionDomesticationResult
TranscriptChoice
DesignConfiguration
```

The fixed Addgene #169227 workflow can then become one tested configuration rather than a special case.

---

# 21. Testing priorities

Maintain the existing Tubb5 fixture as a regression test.

Coverage present in the 28-test `0.4.0` suite:

- plus- and minus-strand C-terminal design paths;
- plus-strand N-terminal preview;
- guide PAM destroyed by insertion;
- protospacer retention and 14-nt cutoff behavior;
- blocking mutation required;
- synonymous PAM mutation possible;
- synonymous seed mutation possible;
- no safe synonymous mutation;
- internal coding SapI site;
- internal noncoding/UTR SapI site;
- multiple SapI sites;
- reverse-strand guide blocking plus SapI domestication;
- Ensembl 429, 5xx, connection retry, clean error, and transcript-parent validation;
- exact fixed plasmid assembly, zero residual SapI sites, fusion translation, and GenBank round trip;
- SapI quality-control summaries for resolved and unresolved sites.

Still needed, especially for the next chunk:

1. minus-strand N-terminal preview;
2. no NGG guide within the default window;
3. alternative transcript with a different terminal exon;
4. retired/archive transcript ID;
5. malformed custom `.dna` file;
6. malformed or annotation-poor GenBank/FASTA uploads;
7. custom backbone with incorrect site count or overhang order;
8. duplicate/non-unique assembly overhangs;
9. custom payload with frame error or internal stop;
10. custom payload containing SapI;
11. deliberately residual SapI site in a generalized final plasmid;
12. deliberate final fusion translation mismatch;
13. built-in backbone through generic configuration produces an identical fixed result.

For sequence-critical code, tests should assert exact sequences and coordinates, not only that a result object was returned.

---

# 22. Immediate local development checklist

Completed in the current repository:

- [x] create/use the Python 3.11 `hdr-tag-designer` conda environment;
- [x] run all tests and reproduce the bundled Tubb5 result;
- [x] run live mouse and human C-terminal designs;
- [x] add bounded Ensembl retries, clean errors, and transcript-parent validation;
- [x] generalize coding SapI domestication;
- [x] implement automatic synonymous guide blocking;
- [x] add user-facing SapI site/mutation quality control.

Recommended next order:

1. add `BackboneDefinition`, `PayloadDefinition`, and `AssemblyJunction` models;
2. express Addgene #169227 and its fixed payload using those models;
3. route the existing fixed assembler through the configuration layer and prove exact regression equivalence;
4. add normalized SnapGene/GenBank/FASTA upload parsing with file limits and clean errors;
5. add the custom-backbone configuration and QC panel, initially releasing SapI configurations only;
6. add custom DNA/FASTA payload configuration and junction/translation validation;
7. add custom-backbone/payload failure fixtures from section 21;
8. continue transcript-selection UI, caching, and archive/local-sequence fallbacks;
9. add a validated N-terminal backbone;
10. add locus-specific primer design.

Commands:

```zsh
conda env create -f environment.yml
conda activate hdr-tag-designer
python -m pip install -e .
./run_tests.zsh
streamlit run app.py
```

---

# 23. Important scientific limitations

This application is a computational design aid.

Before experimental use, independently verify:

- transcript choice,
- current genome assembly,
- sample genotype,
- guide sequence,
- PAM,
- guide activity,
- homology-arm sequences,
- all introduced mutations,
- cloning extensions,
- type-IIS digest behavior,
- donor-plasmid assembly,
- 5′ and 3′ edited junctions,
- fusion reading frame,
- tag orientation,
- biological suitability of the selected terminus.

A sequence-valid design is not automatically biologically valid. Terminal tagging can disrupt localization signals, interaction motifs, post-translational modification sites, degradation signals, and protein function.

---

## 24. Reference material used for the prototype

- Bollen et al. (2022), main paper and supplementary materials describing ITPN guide selection, homology-arm design, SapI adapters, and donor architecture.
- Uploaded Addgene plasmid #169227 SnapGene `.dna` file.
- Uploaded Addgene plasmid #169226 SnapGene `.dna` file.
- Uploaded `bollen_supplementary_s1.docx`, used to verify the exact N-terminal arm adapters and primer tails.
- Ensembl REST for live mouse and human transcript/genomic sequence retrieval.
- Bundled mouse Tubb5 reference fixture for regression testing.

---

## 25. Current bottom line

The current version should be considered:

> A sequence-complete N- and C-terminal prototype for reference-based tagging of mouse or human protein-coding transcripts using Bollen/Addgene #169226 or #169227, with a guarded custom SnapGene-backbone path for the same SapI/linker architectures. It automatically applies verified synonymous guide-blocking plus synonymous-coding or guarded non-coding SapI-domestication mutations, exposes those decisions for quality control, and withholds sequence-complete output when no safe automatic solution exists.

The next major extensibility objective is:

1. **split the current generic `BackboneDefinition` into clean serializable backbone and payload definitions;**
2. **add custom DNA/FASTA payload-only input and user tag metadata;**
3. **add GenBank backbone input, comprehensive negative fixtures, and optional circular-origin normalization.**

Preserve these invariants throughout that work:

- the fixed Tubb5 sequence and coordinate regression must remain exact;
- the nearer-guide/silent-mutation preference must not change;
- no order-ready output may be released while any backbone, payload, mutation, junction, frame, topology, or residual-site validation is unresolved.
