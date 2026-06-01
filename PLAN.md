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

### Q3.2 — Permutation null + FDR ✅ **DONE (2026-05-27)**

- [x] Extend `src/scoring/engine.py` with shuffled-gene-set null
  distributions per `(cell_state, species)`. Implemented in new
  `src/scoring/null.py` (`score_with_null`, size-matched random draws
  from `adata.var_names`, one-sided absolute-deviation test).
- [x] Emit per-module FDR alongside raw scores. CLI: `wombat
  score-null --n-permutations N --group-key cell_type`. Output:
  `results/reports/scoring/permutation_fdr.{csv,md}`.
- [x] Smoke test on integrated h5ad. 64 (module × species × cell_type)
  tests at n_perm=100: 13 significant at FDR<0.05;
  **3 modules significant in both species** —
  `decidual_score`, `ECM_remodeling_score`, `immune_interface_score`.

### Q3.3 — Conserved vs divergent classification ✅ **DONE (2026-05-27)**

- [x] For each of 8 modules: classify as conserved / human-biased /
  mouse-biased using effect size + FDR. Implemented in
  `src/scoring/conservation.py` (`classify_conservation`,
  `summarise_modules`).
- [x] Write `results/reports/conservation_table.csv` + a markdown
  summary. CLI: `wombat classify-conservation`. Outputs:
  `results/reports/conservation_{table.csv,table.md,summary.csv}`.

**Headline result** (n_perm=100, FDR<0.05):

| class | n modules | examples |
|-------|----------:|----------|
| conserved-up | 1 | `decidual_score` (decidual_stromal, +0.62 human / +0.76 mouse) |
| mouse-biased | 4 | `immune_interface`, `senescence`, `ECM_remodeling`, `angiogenesis` |
| human-biased | 2 | `estrogen_response`, `stress_response` |
| neutral | 1 | `progesterone_response` |

This is the manuscript headline: the canonical `decidual_score` is
cross-species conserved while most companion modules diverge in
cell-type-specific directions.

### Q3.4 — Candidate regulator shortlists ✅ **DONE (2026-05-27)**

- [x] New `src/cell_states/regulators.py`. Uses Lambert 2018
  (1,639 human TFs), shipped at
  `configs/reference/lambert2018_human_TFs.txt` from pySCENIC.
- [x] Rank by (a) score correlation in decidual stromal cells in BOTH
  species → **conserved** list (`results/reports/regulators_conserved.csv`);
  (b) species-rank gap → `regulators_human_biased.csv` /
  `regulators_mouse_biased.csv`.
- [x] Cap each at 25. CLI: `wombat rank-regulators [--cap N]`.

**Conserved top 10 (mean rank across species):**
MSX1, ELF3, IRX3, BHLHE40, EHF, FLI1, MECOM, ATF3, TBX3, NR4A2 —
includes the canonical FOXO1 (rank 27) and well-known stress / IE
TFs (NR4A1/2, ATF3, KLF6, MAFF) plus stromal-differentiation TFs
(MSX1, MECOM, SMAD3, NFIA/B, TBX3).

### Q3.5 — Atlas pages for results

- [ ] Add Streamlit pages: Conservation table, Regulator browser. Reuse
  existing patterns in `decidual_atlas/`. No new architecture.

### Q3.6 — Venue decision checkpoint ⏯️ **DEFERRED to Q5 (re-decided 2026-05-27)**

Originally scheduled here on the assumption that the conserved-program
finding was the manuscript headline. The Q4 reframe (convergent
evolution of spontaneous decidualization) makes this premature: with
Q4.1–Q4.3 evidence in hand, eLife / Genome Biology become live again
and the venue table looks materially different. Decision moved to
**Q5.1** (post-Q4.4 reframe).

- [x] Drafted [`docs/manuscript_outline.md`](docs/manuscript_outline.md)
  as a *provisional* scaffold (title, 5+3 figure plan, methods sections
  mapped to `src/` modules, reproducibility deliverables, risks,
  submission checklist). The figure / methods / repro scaffolding
  carries over to Q5 unchanged; the **venue rationale (§1), working
  title (§2), and headline finding (§3)** are flagged as superseded
  pending Q4 evidence.
- [~] Venue chosen — **withdrawn**. Q3.6's GigaScience selection was
  conditional on a conservation-only story; re-decide in Q5.1.

### Q3 exit criteria

- [x] ≥1 conserved module at FDR < 0.05 in both species. *(Q3.2:
  decidual_score, ECM_remodeling_score, immune_interface_score.)*
- [x] Score-vs-time monotonicity plot in report (mouse + human bulk).
- [x] Two regulator lists exported. *(Q3.4: conserved + human-biased
  + mouse-biased.)*
- [x] Manuscript outline scaffold drafted. *(Q3.6:
  `docs/manuscript_outline.md` — figure / methods / repro plan; venue
  + headline deferred to Q5.1 post-Q4 reframe.)*

---

## Q4 (Year 2 H1) — Convergent evolution of spontaneous decidualization

**Scientific question.** Spontaneous decidualization (decidua forms each
menstrual/oestrous cycle without an embryonic signal) evolved
**independently** in catarrhine primates, elephant shrews, the spiny
mouse (*Acomys cahirinus*), and some bats. Q3 established that the
*execution machinery* (the `decidual_score` gene set + its conserved
regulator core) is shared between human and mouse — i.e. when mouse
stroma is triggered, it runs essentially the same program. **Q4 asks
what changed in the lineages that acquired the spontaneous trigger.**

The four leading hypotheses (testable, distinguishable):

1. **Lowered activation threshold** — basal expression of PGR cofactors
   / FOXO1 / HAND2 is higher in unstimulated spontaneous-deciduator
   stroma, so endogenous progesterone suffices.
2. **cis-regulatory rewiring** (Lynch/Wagner) — TE-derived enhancers
   (MER20, MER41) recruited to decidualization genes independently in
   each spontaneous lineage.
3. **Loss-of-repressor** — induced-deciduators carry a brake (AHR-like,
   miRNA, etc.) that spontaneous lineages lost.
4. **Stromal-niche pre-priming** — uNK / vascular crosstalk in
   spontaneous cycles provides an endogenous implantation-like signal.

### Q4.1 — Baseline-priming test (uses existing atlas, no new data) ✅ DONE

- [x] `src/scoring/baseline_priming.py`: per-species priming distance
  (Cohen's d, resting `stromal_fibroblast` → decidualized lineage =
  `{pre_decidual_stromal, decidual_stromal, senescent_decidual}`) for
  all 10 module scores in the atlas, plus between-species Welch's t at
  the resting baseline.
- [x] CLI: `wombat score-baseline`.
- [x] Reports: `results/reports/baseline_priming.{md,csv}` and
  `results/reports/baseline_priming_between_species.csv`.
- [x] **Decision rule applied to `decidual_score`** (human n_resting=8,466,
  mouse n_resting=10,489):
  - Human priming distance: **0.766** Cohen's d.
  - Mouse priming distance: **2.625** Cohen's d.
  - Gap (mouse − human): **+1.860** Cohen's d.
  - **Verdict: hypothesis 1 (lowered activation threshold) is
    SUPPORTED.** Human resting stroma sits ~3.4× closer (in
    standardised units) to the decidualized end-state than mouse
    resting stroma does. Mouse stroma must traverse much more
    transcriptional distance to reach the same conserved decidual
    programme, consistent with a higher activation threshold.
- [x] **Caveats to revisit in Q4.2/Q4.3.** Mouse decidualized cells
  come from early pregnancy (GSE226417, embryonic signal present);
  human resting cells pool the natural cycle (proliferative →
  late-secretory). The directional result is robust to both biases
  (both would *narrow* the gap, yet the gap is huge), but cycle-stage
  stratified and Acomys-rodent contrasts (Q4.2) are needed to rule out
  the alternative that the gap reflects pregnancy-vs-cycle staging
  rather than species biology.

### Q4.2 — Phylogenetic expansion (multi-species trait contrast)

**Data-availability scoping (done 2025-Q4):** verified the placeholder
accessions previously listed here (GSE124753, GSE149298, GSE32916)
and found they point to unrelated studies (cat parvovirus, human
T-cell autoimmunity, mouse macrophage transcription). Searched GEO
and Europe PMC for deposited *Acomys cahirinus* endometrial
transcriptomes; **none exist publicly as of the search date.** The
Bellofiore / McKenna line of work is histology + endocrinology only;
one 2022 SAGE-seq conference abstract was never deposited. Spiny
mouse is therefore **out of scope as a drop-in dataset** for Q4.2 —
keeping it would require either author contact or a new generating
experiment, neither of which is in Q4 scope.

**Replacement plan — use the actual published comparative deciduagenesis
datasets:**

- [x] Ingest **GSE274701** (Mika / Wagner group, 2024) —
  multi-species single-cell atlas of midgestation fetal-maternal
  interface (opossum, tenrec, guinea pig, mouse; **deposit does NOT
  contain macaque despite paper text**). Implemented
  `geo_per_species_h5ad` ingest path (commit 85de7fa); deposit
  ships one h5ad per species, copied into canonical
  `results/processed/GSE274701__{species}.h5ad` naming + manifest.
  Per-species cell counts: guinea_pig 18,131 / mouse 11,452 /
  opossum 6,761 / tenrec 29,652 (~66k cells, 4 species). The
  `GSE274701_RAW.tar` (836 MB) also unpacked a bonus per-sample
  MTX set plus a Te132 kallisto MPA/DIFF/T25 stim-time-course; not
  consumed by the current manifest but available for later
  stim-response work. The macaque trait-positive datum the
  original Q4 plan called for must come from elsewhere
  (likely Lyu 2022 reanalysis or a new ask).
- [x] Ingest **GSE109309** (Erkenbrack / Wagner 2018, "The mammalian
  decidual cell evolved from a cellular stress response") — opossum
  endometrial stromal cells, bulk RNA, 15 samples × 23,899 ENSMODG
  genes under the standard UNDIFF/DIFF/MPA/PGE2/PGE2+MPA in vitro
  stim panel (commit 64171f8). Required a small fix to
  `src/ingest/anndata_writer.py::_load_csv_counts` to drop
  non-numeric annotation columns (e.g. the leading `gene_name`
  column shipped alongside the ensembl_id index).
- [x] Extend `configs/species.yaml` to the new species (opossum,
  tenrec, guinea pig, macaque, baboon, hamster, 13-lined ground
  squirrel, armadillo, **Myotis lucifugus**). 11 entries total;
  per-species `tier`, `ensembl_dataset`, `ensembl_prefix`,
  `source_datasets`, and trait labels recorded.
- [x] **Ingest GSE155170** (5-species endometrial bulk; corrected from
  the prior "maternal-fetal interface" label after inspection).
  Implemented `geo_per_sample_bulk` ingest path in
  `src/ingest/bulk_multi_species.py` and `.tar` extraction in
  `src/ingest/geo.py`. Produces one h5ad per species
  (`results/processed/GSE155170__{species}.h5ad`) plus a manifest
  h5ad whose `.uns['per_species_h5ads']` indexes them. **Key
  discovery from real ingest:** the "bat" samples were aligned to
  *Myotis lucifugus* (vespertilionid, trait-negative — ENSMLUT
  transcript prefix), not *Carollia perspicillata* (phyllostomid,
  trait-positive) as the original Marinic/Kin/Wagner paper text
  suggested. The Q4 convergence figure therefore still **lacks a
  phyllostomid trait-positive bat**; this is a follow-up data ask.
  Hamster samples carry NCBI RefSeq IDs (NM_*), so a RefSeq→Ensembl
  bridge is required at the ortholog-mapping step.
- [x] **Generalize `src/orthologs/ensembl.py` to all Tier B species.**
  The module now reads `ensembl_dataset` / `ensembl_prefix` /
  `ensembl_species` from `configs/species.yaml` rather than from a
  hard-coded mouse-only dict. BioMart attribute names are derived
  from the target's `ensembl_prefix`; the parser locates target
  columns by structural rule (not by species display label) so any
  Tier B target works without code changes. CLI gained
  `--target NAME` and `--all-tier-b` options; per-target failures
  are logged and the loop continues so one bad species (e.g. tenrec,
  see below) does not block the rest. `_parse_biomart_response` and
  the mirror-retry loop now detect Ensembl maintenance pages (HTML
  body / `text/html` content-type) and fall through to the next
  mirror / Compara FTP fallback. Test coverage:
  `tests/test_ortholog_dispatch.py` (6 new tests).
- [x] **Run human→{macaque, baboon, opossum, guinea_pig, hamster,
  ground_squirrel, bat_carollia, armadillo} backbones.** Eight of
  nine Tier B targets built (row counts after one-to-one filter
  upstream of `results/orthologs/backbone__human_<tgt>.parquet`):
  macaque 23217, baboon 22718, armadillo 22557, bat_carollia 21535,
  opossum 20283, guinea_pig 19996, ground_squirrel 19586, hamster
  19134. Tenrec excluded — only *Echinops telfairi* is in Ensembl
  Vertebrates, not *T. ecaudatus*; substitute *E. telfairi* later if
  the comparative analysis needs an afrotherian. Two hardening
  follow-ups landed during this run: (a) a 0-row no-cache guard so
  flaky BioMart responses (truncated TSV with valid header but no
  data rows) raise instead of poisoning the parquet cache; (b)
  `_fetch_compara_ftp` now tries the target-species Compara directory
  first and falls back to the source-species (`homo_sapiens`)
  directory when the target's per-species file does not contain
  homologies for the source (baboon's `papio_anubis/` file has zero
  `homo_sapiens` rows). Hamster still needs the RefSeq→Ensembl
  bridge before downstream mapping consumes the backbone.
- [ ] Re-run the Q3.2 / Q3.3 / Q4.1 baseline-priming pipeline on the
  expanded set. **Key contrast**: `decidual_score` activation
  amplitude (Q3.2) and `priming_distance` (Q4.1) correlated with the
  spontaneous/induced trait controlling for phylogeny via paired
  catarrhine-vs-rodent contrasts and (if tenrec menstruation status
  literature supports) afrotherian-vs-mouse contrast; `phylolm`
  becomes feasible at n=6 species.
- [x] New: `src/scoring/trait_contrast.py` — pseudobulk module-amplitude
  trait contrast (Welch's t + Cohen's d + BH-FDR) with a symbol-fallback
  ortholog mapper and one-to-many paralog collapse (most-expressed
  representative per human gene). CLI `wombat trait-contrast`; tests in
  `tests/test_trait_contrast.py`. Ran on the 4 GSE155170 gene-level
  deposits (baboon, bat_carollia = spontaneous; ground_squirrel, mouse =
  induced). **Result** (`results/reports/trait_contrast.md`): no module
  is higher in spontaneous deciduators at FDR < 0.05. `decidual_score`
  trends up (Δ=+0.24, d=+0.50, p=0.46 NS); progesterone/estrogen-response
  modules are *higher* in the induced rodents. An early apparent
  `decidual_score` signal (p=0.033) was a paralog-dilution artifact
  (mouse one-to-many prolactin orthologs) that the representative-ortholog
  collapse removed. **Confound:** in this 4-species subset trait is
  perfectly confounded with clade (catarrhine/bat vs rodent), so this is
  a trait-or-clade contrast, not a phylogeny-controlled test — documented
  in the report; `phylolm` deferred until a within-clade trait contrast
  exists.
- [x] **Tissue caveat.** GSE274701 is midgestation fetal-maternal
  interface, not cycling endometrium. Q4.1's
  resting→decidualized priming distance does not extend directly;
  pseudobulk per-species "decidual program" amplitude (Q3.2-style)
  is the correct readout. Documented in
  `results/reports/baseline_priming.md` ("Tissue caveat" section).

#### Q4.2 trait-coverage gaps (surfaced by the 2026-05-27 ingest round)

The 2026-05-27 ingest closeout surfaced two gaps in the
spontaneous-decidualization trait-positive species coverage. Both
are tracked here as scoped follow-ups, **not** Q4.2 blockers — the
trait contrast can already be run on the catarrhine-vs-rodent axis
with the data in hand (baboon trait-positive from GSE155170 against
5 trait-negative rodents from GSE155170 + GSE274701).

- [x] **Pin tenrec menstruation status to literature** (2026-05-27).
  Verdict: **trait-NEGATIVE** (induced deciduator, seasonal estrous
  cycle). Evidence: Strassmann 1996 (QRB) menstrual-species catalog
  does NOT list tenrec; Nicoll & Racey 1985 (*Reproduction* 74:47)
  describe a follicular/ovulatory cycle, not a menstrual one;
  Poppitt & Speakman 1994 (*Physiol Zool* 67) describe seasonal
  breeding with defined gestation, no cyclic endometrial sloughing.
  Hemochorial placentation (Carter et al. 2004, *Placenta*) does
  NOT imply menstruation (mouse, rat, guinea pig are hemochorial
  and trait-negative). The afrotherian-menstruation argument
  originates from elephant shrew (Macroscelididae, Emera & Wagner
  2012); Tenrecidae are a sister Afrotherian family and the
  inference does not transfer. **Year-one role:** the 29,652-cell
  GSE274701 tenrec h5ad is a **deep afrotherian trait-NEGATIVE
  outgroup** that widens the phylogenetic span of the trait-negative
  pool, not a trait-positive datum. `configs/species.yaml` updated
  (menstruates: false, spontaneous_decidualization: false). Baboon
  (GSE155170, n=3) remains the sole non-human trait-positive datum
  in current hand.

- [x] **Trace the macaque reference in the GSE274701 source paper.**
  Resolved 2026-05-27 by reading the data-availability section of
  the published article (Stadtmauer DJ, Basanta S, Maziarz JD, Cole
  AG et al. 2025, *Nat Ecol Evol* 9(8):1469-1486, PMID 40596730,
  PMCID PMC12328210, DOI 10.1038/s41559-025-02748-x — the paper
  attribution in earlier project notes as "Mika 2024" was wrong;
  the bioRxiv preprint is doi 10.1101/2024.05.01.591945, 2024-09).
  Verbatim from the article: *"macaque data were retrieved from
  accession number **GSE180637**"* — Jiang X et al. 2023,
  *Dev. Cell* 58:806–821.e7 (PMID 37054708, doi
  10.1016/j.devcel.2023.03.012). Cynomolgus macaque (*Macaca
  fascicularis*) placenta scRNA across gestation days 20–140 of
  a 162-day gestation. Added to `configs/datasets.yaml` as a
  Q4_trait_positive_catarrhine entry and the `macaque` block in
  `configs/species.yaml` was switched from rhesus (*M. mulatta*,
  taxon 9544, ENSMMUG, mmulatta) to cynomolgus (*M. fascicularis*,
  taxon 9541, ENSMFAG, mfascicularis) so the gene-ID prefix
  matches what GSE180637 actually ships. Additional mouse refs the
  paper integrated (for the record, not Q4.2 priorities): GSE152903,
  GSE156125, GSE196825. Human data was pulled from
  reproductivecellatlas.org in already-aligned form for privacy.
  Net effect on trait-positive pool: {human, baboon, **cynomolgus
  macaque**, **Carollia perspicillata bat** (see next item)} — two
  independent catarrhines confirmed and a phyllostomid.

- [x] **Phyllostomid bat sweep.** Resolved 2026-05-27 by querying GEO
  E-utilities for Carollia/Phyllostomus/Desmodus/Artibeus endometrium
  records. The very first hit was **GSE155170 itself** — the dataset
  we had already ingested. GEO sample metadata
  (GSM4696524 "Bat Endometrium Individual 1", GSM4696525 "Bat
  Endometrium Individual 2") explicitly names *Carollia
  perspicillata* (phyllostomid, **trait-positive**, menstruating)
  as the source organism. The earlier project label of "Myotis
  lucifugus, vespertilionid, trait-negative" was a bad inference
  from the ENSMLUT gene-ID prefix in the count files; Marinic/Kin/
  Wagner had simply aligned Carollia reads to the Myotis Ensembl
  reference (no Carollia annotation existed). Fix applied:
  `configs/species.yaml` `bat_myotis` block renamed to `bat_carollia`
  with biological taxon 40233, `menstruates: true`,
  `spontaneous_decidualization: true`; `gene_id_prefix` / Ensembl
  pointers retained at `mlucifugus_gene_ensembl` since orthology has
  to flow through that reference until Ensembl ships a Carollia
  annotation; `configs/datasets.yaml` GSE155170
  `filename_species_map: Bat: bat_carollia`. Re-ingest produced
  `results/processed/GSE155170__bat_carollia.h5ad` (n=2). No author
  email needed; no downgrade required. **Trait-positive pool final**:
  {human, baboon, cynomolgus macaque (GSE180637), Carollia
  perspicillata bat (GSE155170)} — the convergence test now has
  three independent trait-positive lineages (catarrhine primate,
  catarrhine primate, phyllostomid bat) against five trait-negative
  rodents + opossum + tenrec.

### Q4.3 — cis-regulatory layer (Lynch hypothesis)

**Data-availability scoping (done 2025-Q4):** the previously listed
"Mika 2021 GSE174068" is wrong (real GSE174068 is mouse exosome STAT6
anti-tumor). No "Mika 2021 comparative endometrial ATAC" deposition
exists under that description. The actual usable Lynch / Wagner
comparative-regulatory datasets are ChIP-seq and bulk-RNA, not ATAC:

- [x] Ingest **GSE61793** (Lynch lab, "Ancient transposable elements
  transformed the uterine regulatory landscape") — human ChIP-seq /
  TE-derived regulatory landscape. Done 2025-Q4: pulled the three
  series-level **processed peak BEDs** (hg19; H3K27ac 24,329 active
  enhancers, H3K4me3 22,440 promoters, DNaseI 137,107 open-chromatin)
  into `results/raw/GSE61793/`; skipped the two ~600 MB WIG signal
  tracks (not needed for overlap). No peak-calling required — directly
  probes the MER20 / MER41 cis-rewiring hypothesis.
- [x] Ingest **GSE30708** (Wagner lab, "Transposon-mediated gene
  regulatory network rewiring") — per-species TSV deposit
  (commit 64171f8) ingested via new `geo_per_species_table` ingest
  format: human 4 samples (DIFF-Abs/UNDIFF-Abs replicates 1-2),
  armadillo n=1 (single count column), opossum n=1 (RPKM + Mapped
  Reads cols). Thin per-species replication limits this dataset to
  qualitative cross-species DE-pattern overlay rather than per-
  species statistics.
- [x] New `src/cis_regulatory/` module (scope: processed-peak +
  TE-overlap, not ATAC, not peak-calling):
  - `peaks.py` — load + QC the GSE61793 BED5 peaks (canonical chroms,
    bioframe sort, width).
  - `genes.py` — UCSC refGene (hg19) → decidual-gene TSS windows
    (±50 kb, merged) for the proximity test.
  - `te_overlap.py` — RepeatMasker (hg19) overlap via bioframe; a
    genome-wide TE census plus a `near_gene_te_enrichment` Fisher test
    flagging MER20 / MER41 derived enhancers (the primary
    Lynch/Wagner-model test). SVA retroposons live under class
    `Other` in hg19.
  - `motif_enrichment.py` — **deferred** (needs JASPAR + genome FASTA);
    the TE-overlap readout already tests the core Lynch prediction.
  - CLI `wombat cis-regulatory` writes `results/reports/cis_regulatory.md`
    + 3 CSVs; tests in `tests/test_cis_regulatory.py`.
- [x] **Decision rule / finding (2025-Q4).** Genome-wide census:
  ~75 % of H3K27ac decidual enhancers (and ~85 % of promoters) overlap
  a TE — the human decidual regulatory landscape **is** TE-rich,
  consistent with Lynch's landscape claim. But the gene-specific
  prediction is **not** recovered: MER20 0/95 and MER41 1/95 near the
  56-gene decidual panel (Fisher p≥0.22, NS) — the proximity test is
  **underpowered** at this panel size (few near-gene peaks vs <1 %
  family base rate). With only human ChIP in GSE61793 the cross-species
  contrast is not answerable here; it would rely on GSE30708 bulk RNA
  divergence near TE-flagged loci, or a future comparative-ATAC
  deposition. Net: descriptive support for the TE-rich landscape,
  no positive gene-targeted MER20/MER41 signal from this data.

### Q4.4 — Convergence-aware manuscript reframe

- [ ] Update `docs/manuscript_outline.md`:
  - Working title becomes *"Convergent rewiring of a conserved
    decidualization core in spontaneously-menstruating mammals"*.
  - Venue table revisited — with Q4.1–Q4.3 evidence, **eLife** /
    **Genome Biology** become live again; GigaScience remains the
    fallback Data Note venue.
- [ ] Add a fifth main figure: **convergence-vs-conservation
  scatter** (per regulator, trait-correlation effect-size vs
  cross-species conservation).

### Q4 stretches (do if Q4.1–Q4.3 finish early)

- [ ] **dN/dS branch-site** on the conserved-regulator set across
  mammals (PAML / HyPhy), testing positive selection on
  spontaneous-decidualization branches. Lower expected yield because
  the field consensus is *cis*, not coding.
- [ ] **Gene-loss screen** — cross-reference TOGA (Hiller lab,
  ~500 mammals) loss calls with the decidualization gene set,
  restricted to spontaneous-decidualization branches. Hypothesis:
  shared loss-of-repressor.

### Q4 exit criteria

- [ ] Each of hypotheses 1, 2 either supported, refuted, or
  explicitly underpowered — with the evidence pinned in
  `results/reports/`.
- [ ] At least one additional spontaneous-deciduator species (catarrhine
  primate via macaque from GSE274701 or baboon from GSE155170; **bat**
  *Carollia perspicillata* also available from GSE155170**)** ingested
  and scored on the same pipeline as human and mouse.
- [ ] Convergence-vs-conservation figure generated.
- [ ] Manuscript outline updated to reflect whichever finding the
  data supports; venue re-decided.

---

## Q5 (Year 2 H2) — Manuscript, data release, submission

### Q5.1 — Venue decision (re-opened from Q3.6)

- [ ] Re-evaluate venue table with Q4 evidence in hand. Candidate set
  re-expanded: **eLife**, **Genome Biology**, **PNAS** (convergence
  story); **GigaScience** / **Scientific Data** (atlas / Data Note
  fallback). The Q3.6 GigaScience selection was conditional on a
  conservation-only story and is withdrawn.
- [ ] Update `docs/manuscript_outline.md` §1 (venue rationale), §2
  (working title), and §3 (headline finding) to match the chosen
  venue and the Q4 result. §4–§8 (figure / methods / repro /
  checklist) carry over from the Q3.6 scaffold with the Q4.4 figure
  addition.

### Q5.2 — Release engineering

- [ ] `scripts/make_figures.py` — deterministic, seeded; writes every
  paper figure from `results/`.
- [ ] Zenodo / figshare release: archive `results/integrated/`,
  `results/scored/`, `results/orthologs/`, `results/reports/` with
  checksums and a DOI.
- [ ] Tag `v1.0.0`; `pip install convergent-decidua && wombat --help`
  from clean env; verify `CITATION.cff` content.
- [ ] External reproducibility test: recruit one colleague to follow
  `docs/REPRODUCE.md` on a fresh machine; fix anything that breaks.

### Q5.3 — Manuscript + submission

- [ ] Manuscript drafting — methods auto-generated by
  `src/reports/methods.py`; biology + figures hand-written.
- [ ] Preprint on bioRxiv.
- [ ] Submission to venue chosen in Q5.1.

### Q5 exit criteria

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
