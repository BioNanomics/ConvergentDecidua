# Manuscript outline — ConvergentDecidua

**Status:** outline drafted 2026-05-27 (Q3.6 venue checkpoint).
**Chosen venue:** **GigaScience** (Data Note / Research Article hybrid).

---

## 1. Venue decision

| Venue | Fit | Decision |
|-------|-----|----------|
| **GigaScience** | Atlas + reusable methods + biological finding all welcome. Allows a Data Note that ships both the harmonised h5ad and the full analysis stack (`wombat` CLI, `src/scoring/{engine,null,conservation}.py`, regulator pipeline). | **Selected.** |
| Scientific Data | Pure data-resource. Would force the conservation/null framework into a separate methods paper, doubling the publication overhead. | Rejected. |
| Genome Biology | Higher impact but tighter biological-novelty bar; our headline ("only decidual_score is cross-species conserved") may read as expected rather than surprising. | Rejected. |
| eLife | Comparable to Genome Biology; the cross-species atlas alone isn't enough without functional follow-up. | Rejected. |

**Why GigaScience**

1. The deliverable is **both** a resource (integrated h5ad, 24,727
   cells × 11,507 genes, two species) and a method (permutation-null
   FDR + conserved/divergent classification + regulator ranking
   reusable on any module + species set).
2. The venue explicitly supports tooling + reproducibility; the
   `wombat` CLI maps cleanly to their "research object" expectations.
3. Citable DOIs for both the atlas and the codebase (Zenodo) fit
   their submission flow.

---

## 2. Working title

> *"A cross-species single-cell atlas of endometrial decidualization
> reveals a conserved transcriptional core embedded in cell-type-
> specific divergent machinery"*

---

## 3. Headline finding

Among 8 hypothesis-driven gene-set modules (`decidual_score`,
`progesterone_response_score`, `estrogen_response_score`,
`stress_response_score`, `senescence_score`, `immune_interface_score`,
`ECM_remodeling_score`, `angiogenesis_score`), **only the canonical
`decidual_score` is conserved at FDR < 0.05 in both human and mouse
decidual stromal cells**. Every companion module diverges in a
cell-type-specific direction (see `results/reports/conservation_table.md`).

The conserved-regulator shortlist recovers the canonical FOXO1 axis
plus stress / immediate-early (NR4A1/2, ATF3, KLF6, MAFF) and
stromal-differentiation (MSX1, MECOM, SMAD3, NFIA/B, TBX3) TFs;
the divergent shortlists identify species-restricted candidate
regulators (~25 each direction).

---

## 4. Figure plan (target: 5 main + 3 supp)

- **Fig 1 — Atlas overview.** Cohort table (4 GSEs covering
  human cycle + decidua + mouse in-vivo + mouse bulk in-vitro time
  course), UMAP (illustrative; LISI ≈ 1 acknowledged), per-species
  cell-type breakdown. Pull data from `results/reports/manifest.md`
  and `results/integrated/stromal_cross_species.h5ad`.
- **Fig 2 — Bulk validation.** Score-vs-time line plots for
  GSE226429 (mouse, Control→Day5, ρ=+1.00) and GSE104721 (human,
  Day0 vs Day4, ρ=+0.88). From
  `results/reports/scoring/bulk_scoring_report.md`.
- **Fig 3 — Conservation heatmap.** Module × (species × cell_type)
  effect-size heatmap with FDR-significance asterisks. From
  `results/reports/scoring/permutation_fdr.csv`.
- **Fig 4 — Conservation summary.** Bar/diverging chart of the 8
  modules tagged conserved / human-biased / mouse-biased / neutral.
  From `results/reports/conservation_summary.csv`.
- **Fig 5 — Regulator landscape.** Scatter of human_rank vs
  mouse_rank for the 1,639 Lambert-2018 TFs, with the three
  shortlists highlighted; companion mini-panels for the top
  conserved (FOXO1, MSX1, MECOM, ATF3, …) showing per-species
  Spearman correlation with `decidual_score`.
- **Supp 1** — Pre-Q3 acceptance gate: integration QC,
  protected-core marker recovery, IGFBP1 diagnostic.
- **Supp 2** — Module composition and species-overrides
  (`configs/markers.yaml`), ortholog backbone Tier 1 stats.
- **Supp 3** — Full conservation table + full regulator ranking.

---

## 5. Methods sections

1. **Datasets and harmonisation** — registry, ingest scripts
   (`src/ingest/`), ortholog backbone (`src/orthologs/`,
   `results/orthologs/backbone.parquet`, Tier 1 1:1 human↔mouse).
2. **QC** — per-assay pipelines (`src/qc/{scrna,bulk,scatac}.py`).
3. **Integration** — Harmony on HVGs with full Tier 1 gene space
   retained in `.X` (Gate A geometry/biology split).
4. **Module scoring** — `scanpy.tl.score_genes` wrapped in
   `src/scoring/engine.py`; per-species overrides in `markers.yaml`.
5. **Permutation null + FDR** — `src/scoring/null.py`. Size-matched
   random gene sets, one-sided absolute-deviation empirical p-value,
   BH across the (module × species × group) grid.
6. **Conservation classification** — `src/scoring/conservation.py`.
   Class rules for conserved / divergent / biased / neutral.
7. **Regulator ranking** — `src/cell_states/regulators.py`.
   Lambert-2018 TFs (1,639) ranked by per-species Spearman of
   expression vs `decidual_score` within the decidual lineage.
8. **Bulk validation** — `src/scoring/bulk.py`; monotonicity by
   Spearman vs time-axis parsed from sample labels.

---

## 6. Reproducibility deliverables

- **Code:** GitHub repo + Zenodo DOI for the tagged release.
- **Atlas data:** `results/integrated/stromal_cross_species.h5ad`
  deposited at figshare or Zenodo (gitignored locally).
- **Reports:** `results/reports/*.md` reproduced by the
  `wombat generate-reports` command.
- **CLI:** `wombat` (pip-installable) reruns every step end-to-end
  (`fetch → qc → orthologs → integrate → score-decidua → score-bulk
  → score-null → classify-conservation → rank-regulators →
  generate-reports`).

---

## 7. Risks / caveats to address up-front

- **LISI ≈ 1.00** after Harmony: integration removes batch but not
  species geometry. UMAP is illustrative only; conservation claims
  rest on **matched-state module scores + bulk pseudobulk**, not on
  shared embedding.
- **IGFBP1 0 % in mouse atlas:** documented as real-biology /
  capture-window (early-pregnancy T55–T105 only).
  See `docs/marker_recovery_plan.md`.
- **n_perm=100** for the current snapshot; final figures will use
  n_perm=1000 (Q3 stretch — same code, longer runtime).
- **GSE104721 has only 2 time-points;** Spearman ρ=±1 is therefore a
  *direction* test, not a regression. Acknowledge in figure caption.
- **Bulk pseudobulk-vs-bulk-RNA-seq normalisation:** scored on
  CPM+log1p in both cases (same `qc_bulk` path).

---

## 8. Submission checklist (to track in PR)

- [ ] Tag v0.1.0 → Zenodo DOI.
- [ ] Upload `stromal_cross_species.h5ad` to Zenodo (cite in Fig 1).
- [ ] Generate all 5+3 figures from `results/reports/*` via a
  `scripts/build_figures.py` driver (next milestone).
- [ ] Camera-ready Methods (auto from this outline).
- [ ] CITATION.cff already present.
- [ ] GigaDB-style data descriptor (their template).
