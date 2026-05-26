# Marker recovery plan (pre-Q3 gate item A companion)

**Scope:** the 8 canonical decidualization markers tracked by
`src/reports/integration_qc.py`. Source of truth for the current
state is `results/reports/integration_qc.md`.

## Current state

| Marker  | In joint var | Human % expr | Mouse % expr | Bucket |
|---------|:------------:|-------------:|-------------:|--------|
| PGR     | ❌ | — | — | **protected core** (lost: HVG) |
| FOXO1   | ✅ | 48.2 | 21.8 | **protected core** (OK) |
| HAND2   | ❌ | — | — | **protected core** (lost: HVG) |
| WNT4    | ❌ | — | — | **protected core** (lost: HVG) |
| IGFBP1  | ✅ | 43.3 | **0.0** | **protected core** (var OK; mouse expression problem) |
| IL15    | ✅ | 25.4 | 22.1 | **protected core** (OK) |
| PRL     | ❌ | — | — | exploratory (Tier 2 / paralog) |
| LEFTY2  | ❌ | — | — | exploratory (Tier 2 / paralog) |

## Triage

### Recoverable via gate item A (HVG carveout)

`PGR`, `HAND2`, `WNT4`, `IGFBP1` should be present in the per-dataset
`adata.var_names` after orthology mapping — they are well-known
Tier 1 1:1 orthologs. They fall out at HVG selection, not at the
join step. **Action:** force-include the protected-core panel in
the integrated h5ad's retained feature set regardless of HVG rank.
Implementation lives in `src/cell_states/integrate.py`; the panel
itself lives in a new `protected_core:` block in `configs/markers.yaml`.

### Needs an orthology + expression check (not just HVG carveout)

`IGFBP1` is **in** the joint var set but expressed in 0 % of mouse
cells. Three possibilities, ordered by likelihood:

1. **Mouse symbol-case or alias mismatch.** Verify `Igfbp1` is the
   exact `var_names` symbol post-orthology mapping for GSE226417
   before concluding biology.
2. **Stage mismatch.** GSE226417 captures early pregnancy decidua;
   IGFBP1 induction may be confined to specific pseudobulk windows.
   Confirm with a per-`cycle_stage` pseudobulk before drawing any
   "human-only" inference.
3. **True biological difference.** Only acceptable as a conclusion
   after (1) and (2) are ruled out.

This investigation is **not** a Q2 fix; it is a Q3 prerequisite for
any IGFBP1-based claim. Log it as part of gate item E re-evaluation.

### Exploratory only — do not carry year-one claims

`PRL` has no 1:1 mouse ortholog (mouse expanded the decidual-prolactin
family: `Prl8a2`, `Prl3c1`, `Prl3d1`, `Prl8a1`). `species_overrides`
already handles this for *scoring*, but the cross-species *expression
recovery* story is paralog-ambiguous and will not survive review as a
"conserved" claim. Same risk class for `LEFTY2` (`Lefty1`/`Lefty2`).

**Action:** keep both in the diagnostic report but do not include
either in the year-one protected-core narrative. They become a Q4
follow-up if a paralog-resolution strategy is ready.

## Acceptance for gate item A

- [ ] `configs/markers.yaml` has a `protected_core:` block listing
      `[PGR, FOXO1, HAND2, WNT4, IGFBP1, IL15]`.
- [ ] `src/cell_states/integrate.py` retains the protected-core panel
      in the integrated h5ad regardless of HVG selection (PCA/Harmony
      still run on HVGs only — geometry vs biology separation).
- [ ] Re-run `wombat generate-reports`; the marker-recovery table in
      `results/reports/integration_qc.md` shows `in_joint_var = True`
      for all six protected-core genes.
- [ ] `tests/test_real_data.py::test_protected_core_markers_survive`
      green (regression floor for gate item B).
- [ ] IGFBP1 mouse-expression check (above) logged with a one-line
      finding in this document before closing the gate.
