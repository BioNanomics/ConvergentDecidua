# ConvergentDecidua — Plan (Single Source of Truth)

> **This file is the single source of truth for project planning.** All
> agents (Copilot, future contributors) read and update this file. Do not
> maintain parallel plans in chat memory, session notes, or other files —
> mirror back here when the plan changes.

## TL;DR

One reproducible human–mouse decidualization atlas paper in 12 months.

- **Q1 done (May 2026)**: mouse scRNA flows end-to-end via R/Seurat bridge;
  20,549 joint stromal cells (11,484 mouse + 9,065 human) Harmony-integrated
  and scored. Reproducibility floor in place (honest coverage report,
  sha256 manifest, real_data tests, REPRODUCE.md, Docker image verified).
- **Q2 in progress**: fix mouse marker recall (49% → ≥80%), CI green, joint
  integration upgrade, human scATAC.
- **Q3**: bulk-data score validation, permutation null + FDR, conserved /
  divergent module classification, candidate-regulator shortlists.
- **Q4**: figures, Zenodo DOI, external repro test, preprint, submission.

## Revised priorities for Q2 (based on Q1 findings)

1. **Mouse stromal recall (49%) is the next blocker** — was Q3, promoted
   to Q2.1. Root cause: PRL family + other human markers without 1:1
   mouse orthologs silently drop during human→mouse symbol mapping.
2. **Pre-existing lint debt** breaks the "CI green" reproducibility
   promise. Clean in Q2.2 (~5 errors total).
3. **Python version policy** — Q1 hit 5 `zip(strict=...)` Py3.10+ bugs.
   Decide in Q2 whether to bump `requires-python` to `>=3.10` (recommended)
   or pin to 3.9 properly with a CI matrix.

## Locked decisions (do not relitigate)

- Mouse unblock: R/Seurat bridge for GSE226417 (done).
- scATAC included (human-only, GSE183771).
- Candidate regulators: conserved + divergent lists, both labeled.
- Target venue: flexible; decide month 9 (GigaScience vs Scientific Data).
- **Out of scope**: DeciduaAI, Enformer, scGPT, Geneformer, LINGER, GENIE3,
  bat / spiny mouse, PostgreSQL (DeciduaForge), spatial, perturbation,
  ChIP-seq, FASTQ→counts, cross-species ATAC, novel GRN inference,
  foundation-model fine-tuning. Future thesis chapters.

## Auxiliary workstream — `hormonal_analysis/` (May 2026)

Root-level research workspace at `hormonal_analysis/` for the
comparative species matrix, cycle-hormone curves (human + mouse + rat
core; spiny mouse + primate comparator exploratory), and birth-control
endogenous-impact plot. Deliberately separate from `src/` and
`wombat`; reusable helpers may be promoted later but **must not**
silently expand year-one atlas scope.

- Reproducible core: human, mouse, rat (seeded from open canonical
  sources tagged in `hormonal_analysis/sources.yaml`).
- Contraception Phase 6A (endogenous cycle impact) seeded; Phase 6B
  (drug PK curves) explicitly deferred until DrugBank / FDA-label
  curation is greenlit.
- No CI integration. Validation = `ruff check`, `ruff format --check`,
  and running the three scripts deterministically from a clean
  checkout.

---

## Q1 (months 1–3) — Unblock mouse + reproducibility floor — ✅ DONE

- [x] R/Seurat tooling layer (`docker/Dockerfile`, `scripts/rdata_to_h5ad.R`)
  — Seurat 5 `layer=` API with Seurat 4 `slot=` fallback.
- [x] Ingest GSE226417 via `src/ingest/seurat_rdata.py` shim; route in
  `src/ingest/anndata_writer.py::_load_from_dir`; `configs/datasets.yaml`
  declares `ingest: {format: seurat_rdata, include: [UE_DSC]}`.
- [x] Mouse QC pass — 23,794 → 23,471 cells (clean upstream data).
- [ ] Mouse metadata harmonization + marker overrides — **deferred to Q2.1**
  (this is what causes the 49% stromal recall).
- [x] Fix `src/reports/coverage.py` — derived from `.obs.dataset`, not file
  existence. Verified no more false positives.
- [x] Release-manifest checksums (`src/reports/manifest.py` already had
  sha256/bytes/mtime; verified, no change needed).
- [x] Real-data smoke tests (`tests/test_real_data.py` + `real_data`
  pytest marker, skipped in CI by default).
- [x] `docs/REPRODUCE.md` — external-runnable walkthrough; both local
  (Homebrew R + Seurat) and Docker paths verified.
- [x] Committed and pushed to `main` (commits `b502ef5`, `d773a89`).

### Q1 exit criteria — all met

- [x] `wombat integrate --mode stromal` h5ad contains both species
  (`{human, mouse}`).
- [x] `pytest -m real_data` green (4/4 incl. `test_integrated_includes_mouse`).
- [x] Coverage report row count matches files on disk.
- [x] Every file in `results/` has a sha256 in `manifest.csv`.

### Q1 findings that reshape Q2–Q4

- 🔴 Only ~49% of UE_DSC mouse cells pass the stromal annotation filter
  despite being pre-curated decidual stromal cells. Cause: human markers
  without 1:1 mouse orthologs (PRL family) silently drop.
- 🟡 Python 3.9 is fragile: 5 `zip(strict=...)` bugs surfaced only when
  real data hit the cross-species path.
- 🟡 5 pre-existing lint errors (`UP017`, `SIM108`, `F841`) still red.
- 🟢 R/Seurat bridge architecture works; reusable for any future Seurat
  RData input (E-MTAB-11491 could share it, with caveats).
- 🟢 Docker reproducibility path verified; build context must be repo
  root (`-f docker/Dockerfile .`).

---

## Q2 (months 4–6) — Fix mouse recall, harden floor, joint integration + scATAC

**Goal:** Honest cross-species stromal embedding (mouse recall >80%) plus
human scATAC linked to the same stromal cells. CI green end-to-end.

### Q2.1 — Mouse marker recall (promoted from Q3, top priority) — ✅ DONE

- [x] Add `species_overrides` block to `configs/markers.yaml` for genes
  without 1:1 orthologs (PRL → mouse `Prl8a2/Prl3c1/Prl3d1/Prl8a1`,
  `Lefty1/Lefty2`, `Muc1`, `Klrb1c/Klrb1/Ncr1`).
- [x] Extend `src/cell_states/annotate.py::_prepare_gene_sets` to consume
  `species_overrides` (merge with backbone-mapped symbols before checking
  membership in `adata.var_names`).
- [x] Extend `src/scoring/engine.py` likewise via `set_name` parameter.
- [x] Re-run integration; **mouse stromal yield 11,484 → 15,662 cells
  (48.9% → 66.7% of UE_DSC)**. Below the 80% aspirational target but a
  major improvement. Residual gap is an annotation-strategy issue
  (idxmax over cell-type scores misclassifies 7,470 UE_DSC cells as
  `epithelial_glandular`), tracked as Q2.4 follow-up.
- [x] Add `test_mouse_stromal_recall` to `tests/test_real_data.py` at
  60% floor (regression guard; aspirational 80%).
- [x] Empirically validated: removing the `epithelial_glandular`/
  `epithelial_luminal` `Muc1` adds drops recall back to 49% — `Muc1`
  must stay in the epithelial sets (dilutes over-confident epithelial
  score in stromal cells).

### Q2.2 — Reproducibility floor: CI green (lint debt) — ✅ DONE

- [x] Fix all 11 pre-existing lint errors: 6 `B905` (zip(strict=)),
  3 `UP017` (`datetime.UTC`), 1 `SIM108` ternary, 1 `F841` unused-var
  (manually fixed in `src/orthologs/ensembl.py` — the ruff `--unsafe-fixes`
  autofix produced a buggy no-op `RuntimeError(...)` constructor call;
  replaced with proper `logger.warning` and removed the dead variable).
- [x] Bump `requires-python = ">=3.11,<3.13"` (was `>=3.9`). Rationale:
  matches `ruff target-version = "py311"` and CI's `python-version: "3.11"`.
  Eliminates the `zip(strict=...)` Py3.9 compat bug class permanently.
  Local venv on 3.9 must be recreated for installs; ruff still runs.
- [x] `ruff check .`, `ruff format --check .`, `pytest`, `pytest -m real_data`,
  and `wombat validate-config` all green locally.
- [x] CI workflow `.github/workflows/ci.yml` already runs all three jobs
  (lint, test, validate-configs) on Python 3.11. Verified.

### Q2.3 — Joint integration upgrade — ✅ DONE

- [x] Add `--orthology-tier {1,12}` flag to
  `src/cell_states/integrate.py` (default 1) and plumb through the
  `wombat integrate` CLI. `12` includes Tier 2 orthogroups; `1` keeps
  the conservative 1:1 default.
- [x] Switch `batch_key` from `['species']` to `['species','dataset']`
  when both vary (multi-species AND >1 dataset per species). Falls
  back to `'species'` or `'dataset'` for the simpler cases.
- [x] Save canonical output as
  `results/integrated/stromal_cross_species.h5ad`; legacy
  `stromal_harmony.h5ad` is now a symlink to it (copy fallback on
  filesystems without symlink support) for back-compat with existing
  scripts and the real_data tests.
- [ ] scVI path (`--method scvi`) — deferred. Path exists, but scvi-tools
  pulls in torch and is GPU-friendly only. Exercise once when GPU is
  available; not a Q2 blocker.
- [x] Verified on real data: 24,727 joint cells, Harmony converges in
  2 iterations, batch variables = `['species', 'dataset']`,
  real_data suite 5/5.

### Q2.4 — Integration QC + annotation rework — ✅ DONE

- [x] **Hierarchical lineage annotation.** Added a `cell_type_lineages`
  block to `configs/markers.yaml` and a two-pass assignment in
  `src/cell_states/annotate.py`: each cell first gets a `lineage`
  (stromal / epithelial / perivascular / endothelial / immune / other)
  by max-over-constituent-cell-type-scores, then the fine-grained
  `cell_type` is picked within the winning lineage. Avoids the
  vote-splitting failure mode where four stromal sub-types lose to a
  single epithelial bucket.
- [x] **Honest recall target.** GSE226417 is "Uterine **Epithelial**
  AND Decidual Stromal" by design. The ~33% non-stromal fraction is
  ~7,500 cells that score genuinely epithelial (Muc1+Krt18+Epcam) —
  not annotation failure. Mouse recall holds at 66.7%
  (15,662 / 23,471) and the regression floor stays at 60%.
- [x] **`src/reports/integration_qc.py`** — new module that emits
  `results/reports/integration_qc.md` with: (a) LISI mixing on
  `obsm['X_pca_harmony']` for `species` and `dataset`, with a
  qualitative mixing label, (b) per-dataset / per-lineage /
  per-cell_type composition table, (c) canonical-marker recovery
  (in-joint-var flag + per-species fraction expressing). Uses
  `harmonypy.lisi.compute_lisi`; subsamples to 5k cells for speed.
- [x] **Wired into `wombat generate-reports`.** Added `tabulate` to
  the core dependency list (required by `pandas.DataFrame.to_markdown`).
- [x] **Findings surfaced** (the report did its job — see new Q3
  risks below):
  - LISI ≈ 1.00 on both `species` and `dataset` → effectively zero
    cross-species mixing in the current Harmony embedding.
  - 5 / 8 canonical decidual markers (PGR, HAND2, WNT4, PRL, LEFTY2)
    are missing from the joint var set — dropped during HVG selection.
  - Follow-ups tracked as Q3 risks, not Q2 blockers.

### Q2.5 — scATAC (GSE183771, human-only) — 🟡 PARTIAL (scope cap activated)

- [x] **Sparse TF-IDF** — `src/qc/scatac.py::_tfidf` rewritten to
  stay sparse end-to-end. The old version called `X.toarray()` on the
  full peak matrix, which OOMs on real scATAC (~10^4 cells ×
  ~10^5 peaks ≈ 20 GB dense). Now uses `scipy.sparse.diags` for
  row/column scaling. Locked in with unit tests.
- [x] **Gene-activity matrix** — new `src/qc/scatac.py::gene_activity`
  function. Signac-style peak-to-gene aggregation: for each gene,
  sum counts of all peaks whose midpoint falls within `[TSS-upstream,
  TSS+downstream]` (window flipped for minus-strand genes). Pure
  scipy/pandas, no Signac dep. Locked in with unit tests covering
  window logic, strand handling, and input-validation errors.
- [x] **Unit tests** in `tests/test_scatac.py` (4 tests, all on tiny
  synthetic matrices — no real data needed).
- [ ] **Fetch + process GSE183771** — DEFERRED to Q3 stretch per the
  scope cap in this section. Downloading + processing fragment-level
  scATAC (~10s of GB raw → cell × peak matrix → QC h5ad) is multi-day
  work; better done with the infrastructure already validated. The
  ingest/QC/integration code paths are ready when the data lands.
- [ ] **Co-embed with stromal RNA** — also deferred to Q3 (depends on
  fetched data).
- [x] **Scope cap activated:** scATAC moved to Q3 stretch as the
  Q2 plan explicitly permitted.
  Do NOT extend Q2.

### Q2.6 — Snakemake DAG cleanup — ✅ DONE

- [x] Rewrote `workflows/Snakefile` to load `configs/datasets.yaml`
  manually (it is a top-level YAML **list**, which Snakemake's
  `configfile:` directive does not accept — that broke `snakemake -n`
  outright). Now exposes `DATASETS`, `SCRNA_ACCESSIONS`,
  `ALL_ACCESSIONS`, and `SPECIES_OF` as module-level constants the
  rule files consume.
- [x] Updated rule files: `qc.smk` uses `SPECIES_OF[acc]` instead of
  re-scanning `config`; `integrate.smk` uses `SCRNA_ACCESSIONS` and
  writes the Q2.3 canonical `stromal_cross_species.h5ad` (not the
  legacy `stromal_harmony.h5ad`); `reports.smk` references the
  integrated h5ad and declares all six report outputs (was missing
  qc_summary / orthologs / integration_qc).
- [x] `rule all` now targets the canonical Q2 artifacts (integrated
  h5ad + all reports) rather than the stale `results/registry.parquet`
  it pointed at before.
- [x] Added a `validate-workflow` job to `.github/workflows/ci.yml`
  that runs `snakemake -n --snakefile workflows/Snakefile --forceall`
  on every push. Catches broken `include:` paths, rule syntax errors,
  and dangling references **without** running any data pipelines.
- [x] Documented the DAG in `docs/dag.md` (Mermaid + per-rule artifact
  table + how to regenerate `docs/dag.dot`). The DOT file is checked
  in for reviewers without graphviz; rendering to SVG is one
  `dot -Tsvg` command.
- [x] Verified locally: `snakemake -n --forceall` reports a 12-job
  DAG (4× fetch + 4× qc + 1× orthologs + 1× integrate + 1× reports
  + all), no errors.

### Q2.7 — Optional: more mouse cells — ⏭️ SKIPPED

- Mouse stromal cell count after Q2.1+Q2.4 is **15,662**, well above
  the 8K threshold that would have triggered this section.
- UE_EC subset and E-MTAB-11491 remain available for Q3 if needed.

### Q2 exit criteria

- [x] One joint h5ad with mouse stromal recall at 66.7% (15,662 /
  23,471 UE_DSC cells). The original 80% target was based on the
  incorrect assumption that GSE226417 is pure stromal; it is
  "Uterine **Epithelial** AND Decidual Stromal" by design — see
  Q2.4 for the recalibration.
- [x] `results/reports/integration_qc.md` produced. LISI mixing is
  currently ≈ 1.00 (zero cross-species mixing) — surfaced honestly by
  the report and tracked as a Q3 risk (theta_species bump + reserved-
  marker HVG carveout).
- [x] FOXO1 / IGFBP1 / IL15 recovered in the joint var set across
  both species. PGR / HAND2 / WNT4 / PRL / LEFTY2 lost to HVG
  selection — Q3 fix.
- [x] Human scATAC gene-activity **infrastructure** ready
  (`src/qc/scatac.py::gene_activity`, sparse-safe TF-IDF, unit tests).
  GSE183771 **fetch + processing** deferred to Q3 stretch per the
  Q2.5 scope cap.
- [x] `ruff check .` green locally; `validate-workflow` + `lint` +
  `test` CI jobs green on push (verified on `530e9d7`).
- [x] `pytest -m real_data` 5/5 (original 4 + `test_mouse_stromal_recall`).
  `pytest -q` 6/6 (added 4 scATAC unit tests in Q2.5).

---

## Q2 closeout review (May 2026) — honest assessment

Q2 was a substantial tranche but it strengthened **infrastructure,
annotation architecture, and diagnostic rigor** more than it produced
biology. Three structural weaknesses must be addressed before any Q3
conserved/divergent claim can be defended:

1. **Integration geometry is not biologically informative yet.**
   `integration_qc.md` shows species LISI ≈ 1.00 and dataset LISI ≈
   1.00 — explicit "no mixing, separated clusters." Supports an honest-
   diagnostics claim, not a robust cross-species atlas claim.
2. **Joint feature space is too thin for a core-program argument.**
   5 / 8 canonical markers (PGR, HAND2, WNT4, PRL, LEFTY2) are absent
   from the joint var set, dropped by HVG selection. Any "conserved
   core decidualization machinery" statement is currently underpowered.
3. **Orthology layer is weakly externally validated.** `orthologs.md`
   notes g:Profiler confirmed 0 / 16,168 Tier 1 mappings. Backbone is
   not invalidated, but reviewers will not accept "mapped orthologs
   exist → comparative biology is secure" without per-gene spot-checks
   for the genes that carry the narrative.

Additional clarifications:

- **"CI green" means code-quality + workflow-syntax green.** Real-data
  tests are skipped by default in CI (`real_data` marker excluded in
  `pyproject.toml`). End-to-end biological reproduction is not yet
  CI-enforced.
- **The 66.7% mouse recall improvement is an annotation/mapping
  correction, not a biological discovery.** GSE226417 is mixed
  epithelial + stromal by design; the gain comes from `species_overrides`
  + hierarchical lineage gating. Describe as recalibration, not yield.

**Bottom-line:** Q2 establishes a defensible *resource / framework*
trajectory. It does not yet support an *evolutionary-biology* claim.
The safe central message today is: "we built and stress-tested a
comparative atlas framework and surfaced the bottlenecks for rigorous
cross-species decidualization analysis." Saying "we have identified
conserved vs divergent decidual programs" is premature.

## Pre-Q3 acceptance gate (must pass before Q3.1 begins)

Set this gate now so Q3 does not become another tooling pile. Q3.1 work
is **blocked** until all four items are checked.

### Gate item A — Separate geometry from biology in the integrated object

The current pipeline conflates HVG selection (a geometry concern, for
PCA/Harmony) with the gene set retained for downstream biological
interpretation. Fix:

- [x] Build PCA + Harmony on HVGs as today, but **retain the full joined
  Tier 1 gene space** (or at minimum a protected core panel) in the
  integrated h5ad's `.X` / `.raw` so marker recovery, module scoring,
  and pseudobulk comparisons see the full gene space. Touchpoint:
  `src/cell_states/integrate.py`. **Done:** integrated h5ad now
  carries 11,507 genes (full Tier 1 joint space) with HVG geometry
  preserved in `obsm['X_pca'/'X_pca_harmony'/'X_umap']`;
  `uns['hvg_used_for_geometry']` records the 3,003 HVGs.
- [x] Define the **protected core decidual panel** for year-one claims:
  `PGR, FOXO1, HAND2, WNT4, IGFBP1, IL15`. Registered in
  `configs/markers.yaml::protected_core`. Force-included in the
  integrated object regardless of HVG selection.
- [x] Split marker narrative: **core** = the six above (ortholog-clean,
  carry year-one claims). **Exploratory** = PRL family, LEFTY2
  (paralog-expanded / Tier 2; do not let them carry the main claim).
  Documented in `docs/marker_recovery_plan.md`.

### Gate item B — Drop-audit in integration QC report

- [x] Extend `src/reports/integration_qc.py` so the canonical-marker
      table labels each missing gene with the **reason**: `lost_orthology`,
      `lost_inner_join`, `lost_hvg`, or `present`. Done in commit
      after `133a047`. Real-data run confirms all 5 missing canonical
      markers (PGR, HAND2, WNT4, PRL, LEFTY2) are `lost_hvg` — i.e.
      the HVG carveout in gate item A will recover them.
- [x] Add a regression test in `tests/test_real_data.py` that asserts
      all six protected-core markers survive into the integrated analysis
      layer (`adata_integrated.var_names`). Floor, not aspirational.
      Done as `test_protected_core_markers_survive`; skips upstream-
      absent genes (those require an orthology fix, not an integration
      fix) and asserts on recoverable ones only.

- [ ] Write `docs/ortholog_spotcheck.md`: for each gene in the protected
  core panel, list the human Ensembl ID, mouse Ensembl ID, source
  (Tier 1 backbone), and an independent cross-check (Ensembl Compara
  web entry URL + one literature citation confirming 1:1 orthology).
  Two pages max. No code change required; it is a manual evidence
  trail for the manuscript and a reviewer-defense doc.

### Gate item C — Ortholog spot-check memo

- [x] **Orthology-validation portion (substantially closed).**
      `docs/ortholog_spotcheck.md` covers all 6 protected-core genes
      with Ensembl Compara URLs, HGNC + MGI cross-references, and
      landmark mouse + human uterine functional citations. All 6 are
      Tier 1 `ortholog_one2one` (confidence=1). Backbone built
      2026-04-26 via Ensembl Biomart (specific Ensembl release not
      pinned — recorded as a known limitation in the memo).
- [x] **Expression-side residual (IGFBP1 only).** `scripts/diagnose_igfbp1_mouse.py`
      (2026-05-27) confirms `Igfbp1` is **0.03 % expressing in raw
      mouse cells** (pre-orthology) and **0.03 % post-remap** in the
      mouse subset of the integrated h5ad. The pipeline is not losing
      signal — there is essentially no signal to begin with in
      GSE226417's T55–T105 early-pregnancy window. Per-`(orig.ident,
      time)` pseudobulk shows 0–3 reads across every stage (one
      outlier in a 25-cell sample; not defensible as biology). Audit
      tables + verdict appended to `docs/marker_recovery_plan.md`
      under `## IGFBP1 mouse audit`. **Year-one consequence:** IGFBP1
      cannot carry a cross-species claim from GSE226417; Q3 IGFBP1-
      based statements must restrict to the human side or wait for a
      later-pregnancy mouse dataset (E-MTAB-11491, Q3 stretch). The
      `protected_core` panel composition is unchanged so the 0 %
      remains visible as a transparency diagnostic in
      `integration_qc.md`.
- [x] **Per-locus alignment sanity check (automated).** `src/orthologs/synteny.py`
      + `wombat orthologs synteny-check` query the Ensembl REST
      `/homology/symbol/{species}/{symbol}` endpoint for each
      `protected_core` gene against the configured target species and
      record orthology type, % identity, and dN/dS to
      `results/orthologs/synteny_at_core_loci.parquet`. Replaces the
      manual CGV eyeball loop (CGV is pairwise-only and canvas-rendered;
      Ensembl REST is multi-species and structured). Year-one default
      target = `mouse`; spiny-mouse / bat partners can be enabled via
      `--targets mouse,spiny_mouse` once their relevance is gated in
      (see "Out of scope" in the Locked Decisions section above —
      year-one runs human↔mouse only).

### Gate item D — Honest reproducibility statement

- [x] Update `docs/REPRODUCE.md` (and the README "CI" badge section if
      applicable) to distinguish:
  - **Code-quality CI** (lint, format, unit tests, config validation,
    workflow dry-run) — green on every push.
  - **Real-data reproducibility** (`pytest -m real_data`, full
    Snakemake DAG) — runs locally / on a populated `results/`; not
    enforced by CI. List the exact commands a reviewer needs.
  Done in commit `133a047`.

### Gate item E — LISI re-evaluation + Q3 evidence-chain decision

- [x] Re-ran `wombat generate-reports` (2026-05-27) against the
      post-Gate-A integrated h5ad
      (`results/integrated/stromal_cross_species.h5ad`, 24,727 cells
      × 11,507 genes, protected-core preserved). **LISI is still
      ≈ 1.00** on both axes (species: median 1.00 / mean 1.02;
      dataset: median 1.00 / mean 1.02). Geometry-vs-biology split
      did not unlock cross-species mixing in the Harmony embedding.
      All six protected-core markers (PGR, FOXO1, HAND2, WNT4,
      IGFBP1, IL15) now report `in_joint_var = True` in
      `results/reports/integration_qc.md`.
- [x] **Decision (pre-committed in the closeout review, now
      formally adopted):** the integrated UMAP is demoted to an
      illustrative-only figure. **Q3 cross-species evidence chain =
      matched-state module scores (`src/scoring/engine.py`) +
      pseudobulk on the full Tier 1 gene space retained in
      `.X` / `.raw` of the integrated h5ad.** No cross-species claim
      may rest on UMAP neighborhood structure, integration cluster
      identity, or LISI-derived mixing arguments. Module-level FDR
      (Q3.2) and per-stage pseudobulk comparisons become the primary
      evidence carriers.

### Q2 closeout deliverables (companions to the gate)

- [x] `docs/q2_closeout.md` — one page: what Q2 solved, what Q2 only
      diagnosed, what remains blocked. Reuses the bullets from the
      closeout review above; no new analysis. Done in commit `133a047`.
- [x] `docs/marker_recovery_plan.md` — short action note for the five
      missing markers (PGR, HAND2, WNT4, PRL, LEFTY2): which are
      recoverable via Gate item A's HVG carveout vs which need Tier 2
      orthology vs which stay exploratory. Done in commit `133a047`.
      Gate-B drop-audit subsequently confirmed all 5 are `lost_hvg`.

- Stage comparability across datasets (cycle-day matching, pregnancy-
  day matching) before any cross-dataset claim.
- Within-species batch sanity checks before any cross-species claim.
- Trait-positive vs trait-negative species controls before any
  evolutionary claim. (Out of scope for year one; do not pre-commit.)

---

## Q3 (months 7–9) — Scoring validation + conserved/divergent modules

**Goal:** Module-level conserved-vs-divergent calls with FDR, plus the
conservative candidate-regulator shortlists. `species_overrides` moved
earlier, so Q3 starts with statistically clean inputs.

> **Pre-Q3 acceptance gate fully closed (all five items A–E checked,
> 2026-05-27).** Q3.1 may now begin. Residual caveats — both
> *diagnosed and documented*, not blocking the gate:
> 1. **Cross-species mixing is zero in Harmony** (LISI ≈ 1.00 even
>    after geometry/biology split). Mitigated by Gate-E decision:
>    UMAP demoted to illustrative; primary cross-species evidence
>    chain is matched-state module scores + pseudobulk on the full
>    Tier 1 gene space.
> 2. **IGFBP1 is essentially absent (0.03 %) in mouse GSE226417**
>    (orthology + remap confirmed clean; signal is not in the
>    dataset). IGFBP1-based Q3 claims must restrict to human-side
>    or wait for E-MTAB-11491 (Q3 stretch). Other five protected-
>    core markers are unaffected.

### Q3.1 — Bulk-data score validation ✅ **DONE (2026-05-27)**

- [x] Score GSE226429 (mouse bulk in-vitro decidualization time course)
  using same modules; show monotonic `decidual_score` vs. day.
  **Result: ρ=+1.00, p=0 (Control → Day 5).**
- [x] Select one public **human** bulk decidualization series and add to
  `configs/datasets.yaml`. **Selected: GSE104721 (Sato 2018, EOGT
  endometrial stromal cells, Day 0 vs Day 4 with 8-br-cAMP+MPA, six
  siRNA-NT samples).** PLAN.md's earlier suggestions GSE107844 (actually
  aortic dissection) and GSE4888 (Affymetrix, not counts) were retired
  as unsuitable.
- [x] Score human bulk; document monotonicity.
  **Result: ρ=+0.878, p=0.021.**
- See: `results/reports/scoring/bulk_scoring_report.md`,
  `src/scoring/bulk.py`, `scripts/ingest_gse104721.py`,
  commits `dc941c9`, this commit.

### Q3.2 — Permutation null + FDR

- [ ] Extend `src/scoring/engine.py` with shuffled-gene-set null
  distributions per `(cell_state, species)`.
- [ ] Emit per-module FDR alongside raw scores.
- [ ] Smoke test on integrated h5ad.

### Q3.3 — Conserved vs divergent classification

- [ ] For each of 8 modules: classify as conserved / human-biased /
  mouse-biased using effect size + FDR.
- [ ] Write `results/reports/conservation_table.csv` + a markdown summary.

### Q3.4 — Candidate regulator shortlists

- [ ] New `src/cell_states/regulators.py`. Use a curated TF list
  (AnimalTFDB or similar — pick in Q3.4 kickoff, do not invent).
- [ ] Rank by (a) score correlation in decidual stromal cells in BOTH
  species → **conserved** list; (b) human-only signal → **divergent**
  list.
- [ ] Cap each at ~25; document conservative thresholds.

### Q3.5 — Atlas pages for results

- [ ] Add Streamlit pages: Conservation table, Regulator browser. Reuse
  existing patterns in `decidual_atlas/`. No new architecture.

### Q3.6 — Venue decision checkpoint (month 9)

- [ ] Compare artifact strength: data resource (Scientific Data) vs
  atlas+methods (GigaScience). Pick one.
- [ ] Draft `docs/manuscript_outline.md` against chosen venue's
  requirements.

### Q3 exit criteria

- [ ] ≥1 conserved module at FDR < 0.05 in both species.
- [x] Score-vs-time monotonicity plot in report (mouse + human bulk).
- [ ] Two regulator lists exported.
- [ ] Venue chosen, outline drafted.

---

## Q4 (months 10–12) — Manuscript, data release, submission

- [ ] `scripts/make_figures.py` — deterministic, seeded; writes every
  paper figure from `results/`.
- [ ] Zenodo / figshare release: archive `results/integrated/`,
  `results/scored/`, `results/orthologs/`, `results/reports/` with
  checksums and a DOI.
- [ ] Tag `v1.0.0`; `pip install convergent-decidua && wombat --help`
  from clean env; verify `CITATION.cff` content.
- [ ] External reproducibility test: recruit one colleague to follow
  `docs/REPRODUCE.md` on a fresh machine; fix anything that breaks.
- [ ] Manuscript drafting — methods auto-generated by
  `src/reports/methods.py`; biology + figures hand-written.
- [ ] Preprint on bioRxiv.
- [ ] Submission to chosen venue.

### Q4 exit criteria

- [ ] Preprint posted.
- [ ] Submission acknowledged by chosen venue.
- [ ] Zenodo DOI minted.
- [ ] GitHub release tagged.

---

## Risk register (updated after Q1)

| Risk | Status | Mitigation |
|---|---|---|
| Mouse stromal recall (~67% post-Q2.1) | 🟢 Resolved — Q2.4 hierarchical lineage gate + honest target (dataset is mixed UE+DSC by design) | Score-margin / celltypist optional, not blocking |
| Python 3.9 compat bugs | 🟢 Resolved (Q2.2) | `requires-python = ">=3.11,<3.13"` |
| Pre-existing lint debt breaks "CI green" claim | 🟢 Resolved (Q2.2) | All 11 errors fixed; CI now actually green |
| **Cross-species LISI ≈ 1.00** (no mixing in Harmony embedding) | � Resolved by Gate-E decision (2026-05-27): post-Gate-A re-run still shows LISI ≈ 1.00; UMAP demoted to illustrative; Q3 cross-species evidence chain = matched-state module scores + pseudobulk on the preserved full Tier 1 gene space | Re-evaluate only if a later integration method (scVI on GPU, Q4 stretch) materially changes the embedding |
| **Canonical markers dropped by HVG selection** (PGR/HAND2/WNT4/PRL/LEFTY2 absent from joint var set) | � Resolved by pre-Q3 gate item A (geometry/biology split + `protected_core` carveout in `integrate.py`). Protected core (PGR/FOXO1/HAND2/WNT4/IGFBP1/IL15) all present in integrated h5ad. PRL/LEFTY2 remain exploratory. | n/a |
| **IGFBP1 0% expression in mouse** (separate signal flagged by Q2.4 report) | � Diagnosed 2026-05-27 (`scripts/diagnose_igfbp1_mouse.py`): orthology clean (Gate-C synteny), remap clean (0.03 % pre = 0.03 % post), per-`(orig.ident,time)` pseudobulk = capture-floor noise across T55–T105. Dataset-capture / stage-coverage issue, not a pipeline bug | IGFBP1 Q3 claims restricted to human side, or wait for E-MTAB-11491 (Q3 stretch); panel composition unchanged so the 0 % stays visible as a transparency diagnostic |
| **Orthology layer not externally validated** (g:Profiler 0/16168 Tier 1 confirmations) | � Resolved for the protected core by `docs/ortholog_spotcheck.md` (5/6 cleared via Ensembl Compara + HGNC/MGI + functional literature; IGFBP1 ortholog OK but mouse 0 % expression caveat tracked separately) | Broader backbone validation remains a Q4-stretch item; not a year-one blocker |
| **"CI green" overstates reproducibility** (real-data tests skipped in CI) | 🟠 New — surfaced in Q2 closeout review | **Pre-Q3 gate item D:** REPRODUCE.md must split code-quality CI vs real-data repro |
| scATAC slips into Q3 | 🟢 Acceptable | Capped to Q3 stretch; do not extend Q2 |
| Human bulk validation dataset not selected | � Resolved | GSE104721 (Sato 2018) scored 2026-05-27, decidual_score ρ=+0.878 |
| Docker build cache invalidates on R upgrade | 🟢 Low | Multi-stage build if it becomes painful |

---

## MVR 0.1 dataset registry

| Accession | Species | Assay | Status | Role |
|---|---|---|---|---|
| GSE111976 | Human | scRNA-seq | ✅ integrated | Endometrium across natural menstrual cycle |
| GSE127918 | Human | scRNA-seq | ✅ integrated | Decidual pathway / stromal trajectory |
| GSE183771 | Human | scATAC-seq | Q2.5 | Chromatin accessibility across menstrual cycle |
| E-MTAB-11491 | Mouse | scRNA-seq | Q3 stretch | Cycling and decidualizing mouse FRT (644 files) |
| GSE226417 | Mouse | scRNA-seq | ✅ integrated (UE_DSC) | Early pregnancy decidua / uterus |
| GSE226429 | Mouse | bulk RNA-seq | QC'd; scoring in Q3.1 | In vitro decidualization time course |

---

## Appendix A — Module architecture reference

This is the long-standing phase decomposition that maps source modules to
responsibilities. Useful as orientation; the Q1–Q4 plan above is the
operational priority list.

### Module layout

```text
wombat/          # CLI + config loader (Click)
src/             # Analysis modules
  ingest/        # geo, arrayexpress, seurat_rdata, anndata_writer
  metadata/      # harmonize, annotate, audit
  qc/            # scrna, scatac, bulk, pseudobulk
  orthologs/     # ensembl, gprofiler, backbone
  cell_states/   # annotate, subset, integrate
  scoring/       # engine, gene_sets, reports
  reports/       # coverage, manifest, methods, qc_report, ortholog_report
decidual_atlas/  # Streamlit visualization app
configs/         # datasets.yaml, species.yaml, markers.yaml
workflows/       # Snakemake rules (fetch, qc, orthologs, integrate, reports)
scripts/         # standalone helpers (rdata_to_h5ad.R)
```

### Dependency graph

```mermaid
graph LR
    Configs --> Ingest
    Configs --> Orthologs
    Ingest --> Metadata --> QC
    QC --> CellStates
    Orthologs --> CellStates
    CellStates --> Scoring --> Atlas
    Atlas --> Reports
```

### Design notes

- **Harmony vs scVI** — Harmony default (CPU-friendly). scVI via
  `--method scvi` for GPU runs.
- **Data storage** — Large h5ad/parquet in `results/` (gitignored). The
  Snakemake workflow re-derives everything from downloads. Git LFS
  configured in `.gitattributes` for any tracked binary files.
- **Processed matrix availability** — Validate data availability early.
  If a dataset lacks processed matrices, a FASTQ→count alignment step
  is added (not currently needed for MVR 0.1).

### Verification per quarter

`pytest -q` green, `ruff check .` green, `snakemake -n` (dry-run) green,
`wombat validate-config` green.
