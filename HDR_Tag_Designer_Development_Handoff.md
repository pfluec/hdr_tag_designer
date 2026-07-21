# HDR Tag Designer — Development Handoff

**Project version reviewed:** `0.3.1`  
**Handoff date:** 2026-07-20  
**Purpose:** Continue development, debugging, and validation of the local HDR-tagging design tool in a new session or local coding environment.

---

## 1. Project objective

The goal is to build a local Streamlit application that designs CRISPR–Cas9 HDR donor constructs for N- or C-terminal endogenous gene tagging in mouse and human cells.

The initial implementation follows the workflow described by **Bollen et al. (2022)**, especially the supplementary protocol describing:

- guide selection near the intended insertion site,
- homology-arm design,
- placement of guide target sequences outside the homology arms for in-trans paired nicking,
- SapI-compatible Golden Gate assembly,
- sequence validation of the donor and edited locus.

The current prototype is intentionally narrow and uses one fixed, validated donor architecture before adding user-defined backbones and tags.

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
- **Initial fixed donor:** Addgene plasmid **#169227**
- **Initial tag:** C-terminal mNeonGreen
- **Validation gene:** mouse **Tubb5**
- **Validation transcript:** `ENSMUST00000001566`
- **Custom backbones and tags:** add after the fixed-backbone workflow is stable
- **Blocking mutations:** must eventually be designed automatically when the intended edit does not sufficiently disrupt the selected guide target

---

## 3. Current package

The packaged prototype is:

```text
HDR_Tag_Designer_v0.3.1.zip
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
- if the design requires a blocking mutation, final order-ready cloning fragments are currently withheld.

The current code correctly **detects** that a blocking mutation is required, but it does not yet automatically design one for arbitrary genes.

This is one of the highest-priority missing features.

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

Important limitation:

> This correction is hard-coded and validated specifically for the bundled Tubb5 fixture. Generic designs currently detect internal SapI sites but do not automatically generate safe synonymous corrections.

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
- PCR-primer tail templates.

The current primer output only supplies the fixed Bollen cloning tails. It does not design locus-specific PCR annealing regions.

---

## 11. Environment and installation issue

The packaged `run.zsh` currently creates a Python `venv` using whichever `python3` is already available:

```zsh
if [[ ! -d .venv ]]; then
  python3 -m venv .venv
fi
```

The package also currently pins:

```text
streamlit>=1.59,<2
```

This may fail when the system Python is older than the Python version required by the selected Streamlit release.

A dedicated conda environment is preferable.

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

The future `run.zsh` should either:

- assume the conda environment is already active, or
- explicitly create/use a named conda environment.

It should not silently install into an unrelated active Python environment.

---

## 12. Known issue: Ensembl transcript resolution

An observed failure was:

```text
Ensembl REST returned HTTP 500 for /lookup/id/ENSMUST...
```

The app worked when the desired transcript ID was supplied explicitly.

Current practical workaround:

```text
Always provide the exact current Ensembl transcript ID.
```

Likely improvements:

1. add bounded retry logic for Ensembl HTTP 429 and 5xx responses;
2. validate that a supplied transcript belongs to the selected gene;
3. display all current protein-coding transcripts for user selection;
4. use the canonical transcript only as a default, not as a hidden assumption;
5. show transcript version, protein length, CDS length, and terminal exon;
6. add an archive-ID fallback for retired Ensembl IDs;
7. cache successful Ensembl responses locally;
8. give a short user-facing error instead of displaying HTML returned by Ensembl;
9. allow a locally supplied sequence/GenBank record when Ensembl is unavailable.

---

# 13. Priority feature: automatic guide-blocking mutations

## 13.1 Goal

When the intended tag insertion does not sufficiently destroy the selected guide target, the program should alter the donor homology arm so the correctly edited allele is resistant to further nicking.

The mutation must not alter the desired protein or disrupt essential local sequence features.

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

# 14. Priority feature: generic SapI domestication

The current version only automatically fixes the known Tubb5 SapI site.

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
7. verify payload reading frame;
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
- correct reading frame,
- valid linker/tag junction,
- correct stop-codon behavior,
- complete circular-plasmid reconstruction where applicable.

---

# 16. Priority feature: custom tags and payloads

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

## 19.1 PCR primer design

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

Add tests for:

1. plus-strand C-terminal gene;
2. minus-strand C-terminal gene;
3. plus-strand N-terminal gene;
4. minus-strand N-terminal gene;
5. guide PAM destroyed by insertion;
6. protospacer split but PAM retained;
7. blocking mutation required;
8. synonymous PAM mutation possible;
9. synonymous seed mutation possible;
10. no safe synonymous mutation;
11. internal SapI site in coding sequence;
12. internal SapI site in intron/UTR;
13. multiple SapI sites;
14. no NGG guide within the default window;
15. alternative transcript with different terminal exon;
16. Ensembl 429 retry;
17. Ensembl 500 retry and clean error;
18. retired transcript ID;
19. malformed custom `.dna` file;
20. custom backbone with incorrect overhang order;
21. custom payload with frame error;
22. final plasmid containing residual SapI sites;
23. final fusion translation mismatch;
24. GenBank export round trip.

For sequence-critical code, tests should assert exact sequences and coordinates, not only that a result object was returned.

---

# 22. Immediate local development checklist

Recommended order:

1. create the Python 3.11 conda environment;
2. install the project in editable mode;
3. run all existing tests;
4. reproduce the bundled Tubb5 result;
5. run one unrelated mouse C-terminal gene with an explicit transcript;
6. run one human C-terminal gene with an explicit transcript;
7. improve Ensembl retry and transcript-selection behavior;
8. generalize SapI domestication;
9. implement automatic guide-blocking mutations;
10. refactor the fixed backbone into a `BackboneDefinition`;
11. add custom backbone upload;
12. add custom tag/payload configuration;
13. add a validated N-terminal backbone;
14. add locus-specific primer design.

Commands:

```zsh
conda env create -f environment.yml
conda activate hdr-tag-designer
python -m pip install -e .
pytest
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
- Ensembl REST for live mouse and human transcript/genomic sequence retrieval.
- Bundled mouse Tubb5 reference fixture for regression testing.

---

## 25. Current bottom line

The current version should be considered:

> A working prototype for reference-based, C-terminal mNeonGreen tagging of mouse or human protein-coding transcripts using the fixed Bollen/Addgene #169227 architecture, provided the user explicitly selects the correct transcript and manually reviews any design requiring blocking mutations or generic SapI domestication.

The next two most important scientific features are:

1. **automatic guide-blocking mutation design**, and  
2. **generic synonymous SapI domestication**.

The next major extensibility feature is:

3. **user-supplied donor backbones and custom tags/payloads**.
