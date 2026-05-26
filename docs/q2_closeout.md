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
weaknesses that the rest of Q2 did not have scope to repair.

- **Cross-species mixing is effectively zero.** Both
  `species` and `dataset` LISI are 1.00 on the Harmony embedding
  (median over 5 000 cells). "No mixing, separated clusters." See
  `results/reports/integration_qc.md`.
- **5 / 8 canonical markers are absent from the joint var set.**
  `PGR`, `HAND2`, `WNT4`, `PRL`, `LEFTY2` did not survive HVG selection.
  Of the three that did, `IGFBP1` shows **0 %** expression in mouse —
  a separate but related signal of cross-species mismatch.
- **Orthology backbone has no external confirmations.**
  `results/reports/orthologs.md` reports g:Profiler confirmed
  **0 / 16 168** Tier 1 mappings. The backbone is not invalidated, but
  it is not externally validated either.

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
