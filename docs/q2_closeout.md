# Q2 closeout memo

**Period:** Months 4–6 of the 12-month plan.
**Status:** Closed for engineering. Pre-Q3 acceptance gate (see
[`../PLAN.md`](../PLAN.md)) blocks Q3.1 until items A–E pass.

## What Q2 solved

- **Mouse marker recall floor raised.** `species_overrides` in
  `configs/markers.yaml` plus the hierarchical lineage gate in
  `src/cell_states/annotate.py` lifted stromal recall on GSE226417 from
  48.9 % → 66.7 %. A `test_mouse_stromal_recall` regression guard
  protects the gain at a 60 % floor.
- **Reproducibility-floor lint debt cleared.** All 11 pre-existing lint
  errors fixed, `requires-python` bumped to `>=3.11,<3.13`, and CI now
  enforces `ruff check`, `ruff format --check`, `pytest`,
  `wombat validate-config`, and `snakemake -n --forceall`
  (`validate-workflow` job).
- **Joint integration upgraded.** `--orthology-tier {1,12}` flag,
  multi-key Harmony (`['species','dataset']` when both vary), canonical
  output at `results/integrated/stromal_cross_species.h5ad` with a
  `stromal_harmony.h5ad` symlink for back-compat.
- **Annotation strategy rebuilt.** Two-pass lineage → cell-type assignment
  via the new `cell_type_lineages` block in `configs/markers.yaml`.
- **Integration diagnostics shipped.** New
  `src/reports/integration_qc.py` emits
  `results/reports/integration_qc.md` with LISI mixing, per-dataset /
  per-lineage composition, and canonical-marker recovery.
- **scATAC primitives hardened.** `src/qc/scatac.py::_tfidf` rewritten
  to stay sparse end-to-end (the old version OOM'd on real data);
  Signac-style `gene_activity()` added; 4 unit tests on synthetic data.
- **Snakemake DAG repaired.** `workflows/Snakefile` rewritten to load
  the (list-valued) `configs/datasets.yaml` manually; CI now dry-runs
  the 12-job DAG on every push.

## What Q2 only diagnosed (did not fix)

The Q2.4 integration QC report did its job: it surfaced three structural
weaknesses, two of which the pre-Q3 acceptance gate (`PLAN.md`) has
since repaired (gate items A + B) and one of which remains open
(gate item E, embedding mixing).

- ~~**5 / 8 canonical markers absent from joint var.**~~ **Resolved by
  gate item A.** Integrated h5ad now carries the full Tier 1 joint
  gene space (11,507 genes) with HVG geometry preserved in `obsm`.
  All six protected-core markers (PGR, FOXO1, HAND2, WNT4, IGFBP1,
  IL15) recovered with real biological signal: PGR 31.9% human /
  85.5% mouse, HAND2 68.8% / 80.6%, WNT4 30.3% / 60.8%, etc.
  PRL / LEFTY2 remain `lost_hvg` (intentionally — exploratory only,
  paralog-ambiguous).
- **IGFBP1 expressed in 0 % of mouse cells** in the recovered space.
  Separate finding flagged in `marker_recovery_plan.md` for an
  orthology + stage check (likely `Igfbp1` case/alias or stage
  mismatch, not biology) before any IGFBP1-based claim.
- **Cross-species mixing is still effectively zero.** Both `species`
  and `dataset` LISI remain at 1.00 on the Harmony embedding even
  with protected-core in HVGs. Gate item E decision (locked in
  `PLAN.md`): **do not force biology through the integrated UMAP.**
  Q3 evidence chain pivots to matched-state module scores +
  pseudobulk on the preserved full gene space.
- **Orthology backbone is not externally validated**
  (g:Profiler 0 / 16 168 Tier 1 confirmations). Gate item C
  (`docs/ortholog_spotcheck.md`) is the remaining open gate item;
  must land before any comparative-biology claim.

## What remains blocked

- Any "conserved core decidualization machinery" claim (the joint
  feature space is too thin and the embedding has no mixing).
- Any "differentially deployed" or "human-biased vs mouse-biased"
  module claim (same reason; needs Q3.2 FDR on top of a defensible
  embedding or, failing that, a matched-state pseudobulk evidence
  chain — see pre-Q3 gate item E).
- Any evolutionary or mechanistic interpretation (out of year-one
  scope regardless, but explicitly not supported by current evidence).
- End-to-end real-data reproducibility in CI (currently only
  code-quality + workflow-syntax — see `REPRODUCE.md`).

## Honest framing for external communication

The defensible Q2 story is: **"we built and stress-tested a
comparative-decidualization atlas framework and surfaced the
bottlenecks that any cross-species analysis must solve first."**

It is **not**: "we identified conserved vs divergent decidual programs."
That claim requires the pre-Q3 gate to pass first.

## Next action

Execute pre-Q3 gate items A–E in [`../PLAN.md`](../PLAN.md). Do not
open Q3.1 until all five are checked.
