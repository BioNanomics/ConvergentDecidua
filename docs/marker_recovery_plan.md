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

## IGFBP1 mouse audit

Generated 2026-05-27 by `scripts/diagnose_igfbp1_mouse.py`
(Pre-Q3 Gate C residual closeout). Re-run after any re-ingest or
re-QC of GSE226417.

### Q1 — symbol survival (raw mouse symbols, pre-orthology)

| h5ad      |   n_cells | Igfbp1_in_var   | layer   |   pct_expr_all_cells |
|:----------|----------:|:----------------|:--------|---------------------:|
| processed |     23794 | True            | .X      |                 0.03 |
| qc        |     23471 | True            | counts  |                 0.03 |

### Q2 — post-remap survival (integrated joint h5ad)

| view                    |   n_cells | IGFBP1_in_var   |   pct_expr |
|:------------------------|----------:|:----------------|-----------:|
| integrated (all)        |     24727 | True            |      15.91 |
| integrated (mouse only) |     15662 | True            |       0.03 |
| integrated (human only) |      9065 | True            |      43.36 |

### Q3 — per-(sample × time) pseudobulk on QC mouse h5ad

| orig.ident   | time   |   n_cells |   Igfbp1_counts |   lib_size |   Igfbp1_cpm |
|:-------------|:-------|----------:|----------------:|-----------:|-------------:|
| E105         | T105   |        66 |               0 |    1418295 |        0     |
| U105         | T105   |      2572 |               0 |   82320738 |        0     |
| U55.E65      | T55    |      5026 |               3 |  109567761 |        0.027 |
| U65          | T65    |      3231 |               0 |   90422778 |        0     |
| U.E75        | T75    |      3329 |               1 |   84549111 |        0.012 |
| U85.1        | T85    |      1156 |               0 |   43246033 |        0     |
| U85.2        | T85    |      2303 |               2 |   76812098 |        0.026 |
| E95          | T95    |        25 |               1 |     183043 |        5.463 |
| U95          | T95    |      5763 |               1 |   87784509 |        0.011 |

**Verdict: real-biology / dataset-capture.** `Igfbp1` is expressed in
only 0.03 % of mouse cells already at the QC stage (pre-orthology,
mouse symbol). Post-remap mouse-side `IGFBP1` reads 0.03 % — the
remap is *not* losing signal; the signal is not there to begin with.
Pseudobulk by `(orig.ident, time)` shows 0–3 reads per ~3k-cell sample
across every T55–T105 window, i.e. essentially capture-floor noise.
The only outlier (E95: 5.46 CPM) is a 25-cell pseudobulk and not
defensible as biology.

Hypothesis (1) symbol-case is ruled out (`Igfbp1` is present and
identical pre- and post-remap). Hypothesis (3) is consistent with
the data, but cannot be cleanly separated from (2) using GSE226417
alone — the dataset only spans early-pregnancy T55–T105, so a
genuine late-decidualization induction window would be invisible
here regardless. The protected-core *panel* stays intact (we want
this 0 % to remain visible in `integration_qc.md` as a diagnostic),
but **IGFBP1 cannot carry a year-one cross-species claim from
GSE226417.** Q3 IGFBP1-based statements must either (a) restrict
to the human side, or (b) wait for a mouse dataset that covers
later pregnancy / pseudo-pregnancy windows (E-MTAB-11491 is the
queued candidate; flagged in the Q3 stretch lane).

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
