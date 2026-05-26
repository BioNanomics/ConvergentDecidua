# Ortholog spot-check — protected core decidual panel

**Purpose.** Pre-Q3 acceptance gate item C (see [`../PLAN.md`](../PLAN.md)).
The project's ortholog backbone reports g:Profiler confirmed **0 / 16 168**
Tier 1 mappings (see `../results/reports/orthologs.md`). That single-source
posture is not invalidating, but it is not externally validated either. This
memo provides a per-gene, independent evidence trail for the six genes that
underpin year-one comparative-biology claims, so a reviewer can verify each
1:1 human–mouse relationship without re-deriving the backbone.

**Scope.** The six genes in
[`../configs/markers.yaml`](../configs/markers.yaml)::`protected_core`:
`PGR, FOXO1, HAND2, WNT4, IGFBP1, IL15`. PRL and LEFTY2 are intentionally
**out of scope** (paralog-ambiguous; exploratory only — see
[`marker_recovery_plan.md`](marker_recovery_plan.md)).

**What this memo does and does not validate.** This memo validates
*sequence orthology* (1:1 human–mouse mapping with independent
identifier-system cross-references) and *uterine functional relevance*
in both species (landmark loss-of-function or knockdown evidence). It
does **not** by itself validate expression detectability in the
specific datasets the project consumes, stage-matched comparability of
expression across human and mouse, or retention of each gene through
the integration pipeline. Those are separate concerns tracked by
[`marker_recovery_plan.md`](marker_recovery_plan.md), the drop-audit
column in `results/reports/integration_qc.md`, and the
`test_protected_core_markers_survive` regression test.

**Backbone provenance.** Built on **2026-04-26** via Ensembl Biomart
(`src/orthologs/ensembl.py`). The specific Ensembl release in use at
that date is **not pinned** in the project metadata — a known
limitation. If the backbone is rebuilt against a later Ensembl release,
re-run the verification block at the end of this memo and update any
Ensembl gene IDs that have changed. HGNC and MGI cross-references are
stable across Ensembl releases.

**Method.** For each gene we report (a) the human and mouse Ensembl gene IDs
and symbols carried by `results/orthologs/backbone.parquet`, (b) the
Ensembl Compara orthology type and confidence as recorded by the backbone,
(c) the public Ensembl Compara web URL a reviewer can open to inspect the
gene tree and homology relationships, (d) the HGNC ↔ MGI cross-reference
identifiers (an independent identifier system rooted in human and mouse
nomenclature committees, not in Ensembl), and (e) one or two landmark
functional citations confirming the orthologous role in uterine biology.
Where biological evidence is weaker for the mouse counterpart, this is
called out.

**Backbone snapshot** (read from
`results/orthologs/backbone.parquet`, all six rows):

| source_symbol | source_gene_id (human) | target_symbol | target_gene_id (mouse) | orthology_type    | confidence | tier |
|---------------|------------------------|---------------|------------------------|-------------------|-----------:|-----:|
| PGR           | ENSG00000082175        | Pgr           | ENSMUSG00000031870     | ortholog_one2one  | 1          | 1    |
| FOXO1         | ENSG00000150907        | Foxo1         | ENSMUSG00000044167     | ortholog_one2one  | 1          | 1    |
| HAND2         | ENSG00000164107        | Hand2         | ENSMUSG00000038193     | ortholog_one2one  | 1          | 1    |
| WNT4          | ENSG00000162552        | Wnt4          | ENSMUSG00000036856     | ortholog_one2one  | 1          | 1    |
| IGFBP1        | ENSG00000146678        | Igfbp1        | ENSMUSG00000020429    | ortholog_one2one  | 1          | 1    |
| IL15          | ENSG00000164136        | Il15          | ENSMUSG00000031712     | ortholog_one2one  | 1          | 1    |

All six are Ensembl Compara high-confidence `ortholog_one2one` calls. No
Tier 2 / paralog ambiguity in the protected core.

---

## PGR — Progesterone receptor

- **Human:** `PGR`, HGNC:8910, Ensembl `ENSG00000082175` (chr 11q22.1).
- **Mouse:** `Pgr`, MGI:97567, Ensembl `ENSMUSG00000031870` (chr 9).
- **Compara:** one-to-one, confidence 1.
- **Compara URL:**
  <https://www.ensembl.org/Homo_sapiens/Gene/Compara_Ortholog?g=ENSG00000082175>
- **HGNC ↔ MGI cross-ref:**
  <https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:8910> →
  MGI:97567.
- **Functional confirmation in uterus.** *Pgr*-null female mice exhibit
  complete uterine implantation failure and a non-decidualizing stroma
  (Lydon et al., 1995, *Genes & Development* — the foundational
  loss-of-function paper establishing PGR's conserved decidualization
  role). The stromal-specific *Pgr*-cKO phenotype recapitulates the
  decidual defect in mouse (Franco et al., 2012, *FASEB J*),
  supporting functional orthology, not just sequence orthology.
- **Risk:** none. This is the single most uncontested core decidualization
  gene in the panel and safe to centre comparative claims on.

## FOXO1 — Forkhead box O1

- **Human:** `FOXO1`, HGNC:3819, Ensembl `ENSG00000150907` (chr 13q14.11).
- **Mouse:** `Foxo1`, MGI:1890077, Ensembl `ENSMUSG00000044167` (chr 3).
- **Compara:** one-to-one, confidence 1.
- **Compara URL:**
  <https://www.ensembl.org/Homo_sapiens/Gene/Compara_Ortholog?g=ENSG00000150907>
- **HGNC ↔ MGI cross-ref:**
  <https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:3819> →
  MGI:1890077.
- **Functional confirmation in uterus.** FOXO1 is a required transcription
  factor for human endometrial stromal cell decidualization (Christian
  et al., 2002, *JBC* — siRNA knockdown of FOXO1 abolishes IGFBP1 and PRL
  induction in HESC). Uterine-specific *Foxo1* deletion in mouse
  recapitulates a decidualization defect (Vasquez et al., 2018,
  *Endocrinology*), confirming that the decidualization role is conserved
  in the mouse 1:1 ortholog, not just sequence-conserved.
- **Risk:** none for sequence orthology. Cross-species expression-level
  comparison is sound, but **interpretation** of FOXO1 activity in mouse
  must acknowledge that PI3K/AKT regulation of FOXO1 nuclear localization
  may differ in timing relative to human (Vasquez 2018 caveats).

## HAND2 — Heart and neural crest derivatives expressed 2

- **Human:** `HAND2`, HGNC:4808, Ensembl `ENSG00000164107` (chr 4q34.1).
- **Mouse:** `Hand2`, MGI:103580, Ensembl `ENSMUSG00000038193` (chr 8).
- **Compara:** one-to-one, confidence 1.
- **Compara URL:**
  <https://www.ensembl.org/Homo_sapiens/Gene/Compara_Ortholog?g=ENSG00000164107>
- **HGNC ↔ MGI cross-ref:**
  <https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:4808> →
  MGI:103580.
- **Functional confirmation in uterus.** Uterine-specific *Hand2* deletion
  in mouse produces a complete implantation block with sustained epithelial
  FGF signalling and failed stromal differentiation (Li et al., 2011,
  *Science* — the foundational mouse paper). HAND2 is also a validated
  decidualization regulator in primary human endometrial stromal cells
  (Huyen & Bany, 2011, *Reproduction*; Marinić et al., 2021,
  *eLife* — comparative endometrial transcription factor analysis across
  placental mammals explicitly tracks HAND2 as a conserved decidual
  regulator).
- **Risk:** none for sequence orthology; HAND2 is in fact one of the
  textbook examples of a *deeply conserved* decidual transcription factor
  across placental mammals, which strengthens, rather than weakens, its
  use in the protected core.

## WNT4 — Wnt family member 4

- **Human:** `WNT4`, HGNC:12783, Ensembl `ENSG00000162552` (chr 1p36.12).
- **Mouse:** `Wnt4`, MGI:98957, Ensembl `ENSMUSG00000036856` (chr 4).
- **Compara:** one-to-one, confidence 1.
- **Compara URL:**
  <https://www.ensembl.org/Homo_sapiens/Gene/Compara_Ortholog?g=ENSG00000162552>
- **HGNC ↔ MGI cross-ref:**
  <https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:12783> →
  MGI:98957.
- **Functional confirmation in uterus.** Uterine *Wnt4* is required for
  stromal proliferation, differentiation and decidualization in mouse
  (Franco et al., 2011, *FASEB J* — uterine-specific *Wnt4* cKO produces
  a severe decidualization defect). In humans, WNT4 is upregulated during
  the secretory phase and in decidualizing HESC, with WNT4 knockdown
  impairing IGFBP1/PRL induction (multiple primary-cell studies; the
  conserved-decidualization-program review by Gellersen & Brosens, 2014,
  *Endocrine Reviews*, summarises both species).
- **Risk:** low. *WNT4* is part of a larger Wnt-family expression program
  in uterus; co-expression with paralogs (*WNT5A*, *WNT7A*) means that
  module-level interpretation must not over-attribute Wnt-pathway signal
  to WNT4 alone. The 1:1 orthology and per-gene expression call are sound.

## IGFBP1 — Insulin-like growth factor binding protein 1

- **Human:** `IGFBP1`, HGNC:5469, Ensembl `ENSG00000146678` (chr 7p12.3).
- **Mouse:** `Igfbp1`, MGI:96435, Ensembl `ENSMUSG00000020429` (chr 11).
- **Compara:** one-to-one, confidence 1.
- **Compara URL:**
  <https://www.ensembl.org/Homo_sapiens/Gene/Compara_Ortholog?g=ENSG00000146678>
- **HGNC ↔ MGI cross-ref:**
  <https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:5469> →
  MGI:96435.
- **Functional confirmation in uterus.** IGFBP1 is the canonical secreted
  marker of differentiated human decidual stromal cells and is the most
  widely used IVF-era decidualization read-out in human (reviewed in
  Gellersen & Brosens, 2014, *Endocrine Reviews*).
- **⚠️ Cross-species caveat (active investigation — see
  [`marker_recovery_plan.md`](marker_recovery_plan.md)).** The integration
  QC report shows IGFBP1 at 43.4 % expressing in human stromal cells but
  **0 %** in mouse, despite confirmed 1:1 orthology and presence in the
  mouse var set. Three non-exclusive candidate explanations, **not yet
  ranked** — ranking requires the raw-vs-remapped mouse matrix audit
  that has not been performed:
  1. **Symbol / case mismatch survives the backbone remap** in the
     mouse processed h5ad (e.g. raw symbol is `Igfbp1` but our remap
     step is case-sensitive somewhere).
  2. **Stage mismatch.** GSE226417 samples early pregnancy decidua;
     mouse uterine *Igfbp1* induction is reported to be confined to
     specific decidual time points and is much weaker than the
     prolactin-family secretory program in mouse.
  3. **True biological divergence.** The mouse decidualization
     secretome emphasises the *Prl* family (`Prl8a2`, `Prl3c1`,
     `Prl3d1`) rather than IGFBP1 as the dominant secreted product;
     some loss of IGFBP1 dominance in mouse stroma may be biologically
     real, but a strict 0 % value should not be assumed to reflect this
     until (1) and (2) are explicitly ruled out by the audit.
- **Risk:** orthology itself is sound; the 0 % mouse value is a
  pipeline / expression-mapping problem that must be diagnosed before
  IGFBP1 carries any year-one cross-species claim. This is the only
  protected-core gene with an open caveat.

## IL15 — Interleukin 15

- **Human:** `IL15`, HGNC:5977, Ensembl `ENSG00000164136` (chr 4q31.21).
- **Mouse:** `Il15`, MGI:103014, Ensembl `ENSMUSG00000031712` (chr 8).
- **Compara:** one-to-one, confidence 1.
- **Compara URL:**
  <https://www.ensembl.org/Homo_sapiens/Gene/Compara_Ortholog?g=ENSG00000164136>
- **HGNC ↔ MGI cross-ref:**
  <https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:5977> →
  MGI:103014.
- **Functional confirmation in uterus.** IL15 is a conserved decidual
  cytokine essential for uterine NK (uNK) cell recruitment and
  differentiation in both species. *Il15*-null mice lack uNK cells and
  show defective decidual vascular remodelling (Ashkar et al., 2003,
  *J Immunol* and follow-up work by the Croy group). In humans, IL15 is
  produced by decidualizing stromal cells and drives CD56^bright^ uNK
  expansion (Verma et al., 2000, *Hum Reprod*; reviewed in Gaynor &
  Colucci, 2017, *Front Immunol*).
- **Risk:** low for orthology. **Interpretive risk:** IL15 in the
  *stromal* compartment is a signalling source for a *non-stromal*
  effector (uNK); per-cell stromal expression numbers should not be
  read as a proxy for downstream NK-mediated decidual function.
  Module-level use is fine; single-gene cross-species comparison is
  fine; downstream cellular-outcome claims need explicit uNK data.

---

## Summary

| Gene   | 1:1 orthology | Functional support | Present in integrated analysis layer | Open issue | Safe now for year-one claim |
|--------|:------:|:------:|:------:|---|:------:|
| PGR    | ✅ | ✅ (Lydon 1995 mouse cKO + human HESC)            | ✅ (31.9 % human / 85.5 % mouse) | none                                                | ✅ |
| FOXO1  | ✅ | ✅ (Christian 2002 human; Vasquez 2018 mouse cKO)  | ✅ (50.5 % / 43.1 %)             | timing of PI3K/AKT regulation may differ            | ✅ |
| HAND2  | ✅ | ✅ (Li 2011 mouse; Huyen 2011 + Marinć 2021)       | ✅ (68.8 % / 80.6 %)             | none                                                | ✅ |
| WNT4   | ✅ | ✅ (Franco 2011 mouse cKO + human HESC)            | ✅ (30.3 % / 60.8 %)             | family co-expression — module-level only             | ✅ |
| IGFBP1 | ✅ | ✅ in human; mouse role weaker / Prl-family-dominated | ⚠️ present in var, **0 % mouse**  | symbol/stage/divergence audit not yet performed     | ❌ |
| IL15   | ✅ | ✅ (Ashkar 2003 mouse; Verma 2000 human)           | ✅ (25.9 % / 23.8 %)             | downstream cell-fate claims need uNK data           | ✅ |

Integrated-layer percentages are the per-species fraction of cells with
non-zero expression in `results/integrated/stromal_cross_species.h5ad`
after pre-Q3 gate item A's HVG carveout (see
`results/reports/integration_qc.md`). Without gate A, PGR, HAND2, and
WNT4 were `lost_hvg` in the integrated analysis layer despite their
sound orthology — orthology validation alone was never going to make
those markers visible to downstream analysis.

**Bottom line.** Sequence orthology and uterine functional relevance for
the protected core are solid. Of the six, **five (PGR, FOXO1, HAND2,
WNT4, IL15) are safe to centre year-one comparative claims on without
further orthology validation work** — they still must survive whatever
specific analysis (matched-state pseudobulk, module scoring) the
manuscript actually performs, and stage comparability must be
established separately. **IGFBP1** has a 1:1 ortholog and human evidence,
but the mouse 0 % expression needs to be diagnosed (symbol audit →
stage audit → then biological interpretation) before IGFBP1 carries any
cross-species claim.

This memo **substantially discharges the orthology-validation portion of
pre-Q3 acceptance gate item C** in [`../PLAN.md`](../PLAN.md). The
IGFBP1 expression-side audit remains open and is the residual gate-C
work; it is scoped to IGFBP1-based claims only, not a general Q3 blocker.

## How to refresh / verify this memo

```bash
# Re-print backbone rows for the panel
python -c "
import pyarrow.parquet as pq
df = pq.read_table('results/orthologs/backbone.parquet').to_pandas()
core = ['PGR','FOXO1','HAND2','WNT4','IGFBP1','IL15']
print(df[df['source_symbol'].isin(core)].to_string(index=False))
"
```

If Ensembl gene IDs in the backbone change after a future Compara release,
update the Ensembl ID columns and the Compara URLs (the URL format is
stable; only the `?g=` argument changes). HGNC and MGI IDs are designed
to be stable across releases and should not change.
