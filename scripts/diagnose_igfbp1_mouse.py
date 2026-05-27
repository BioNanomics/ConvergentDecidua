"""Diagnose IGFBP1 mouse 0%-expression problem (Pre-Q3 Gate C residual).

Answers three questions, ordered by the hypothesis triage in
``docs/marker_recovery_plan.md``:

1. **Symbol survival** — does ``Igfbp1`` exist in the mouse processed
   and QC'd h5ads, and at what fraction of cells does it have non-zero
   counts?
2. **Post-remap survival** — does ``IGFBP1`` exist in the integrated
   joint h5ad, and at what fraction of the *mouse* cells in the
   integrated object does it have non-zero counts?
3. **Per-stage pseudobulk** — using the QC'd mouse h5ad (the only
   level that still carries the GSE226417 ``time`` column, which is
   the early-pregnancy day proxy), pseudobulk-aggregate ``Igfbp1``
   by ``(orig.ident, time)`` and print per-group summed counts.

The script writes one markdown section to stdout that can be appended
verbatim to ``docs/marker_recovery_plan.md`` (under the existing
"Needs an orthology + expression check" bullet).

Run:

    python scripts/diagnose_igfbp1_mouse.py

No CLI flags. Path-hardcoded to year-one MVR layout; rerun after any
re-ingest or re-QC of GSE226417.
"""

from __future__ import annotations

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

REPO_ROOT = Path(__file__).resolve().parent.parent
PROCESSED = REPO_ROOT / "results" / "processed" / "GSE226417.h5ad"
QC = REPO_ROOT / "results" / "qc" / "GSE226417.h5ad"
INTEGRATED = REPO_ROOT / "results" / "integrated" / "stromal_cross_species.h5ad"

sys.path.insert(0, str(REPO_ROOT))
from src.qc.pseudobulk import pseudobulk  # noqa: E402


def _pct_expressing(adata: ad.AnnData, gene: str, layer: str | None = None) -> float:
    """Return percent of cells with non-zero counts for ``gene``."""
    if gene not in adata.var_names:
        return float("nan")
    j = adata.var_names.get_loc(gene)
    X = adata.layers[layer] if layer and layer in adata.layers else adata.X
    col = X[:, j]
    arr = col.toarray().ravel() if sp.issparse(col) else np.asarray(col).ravel()
    return round(float((arr > 0).mean()) * 100, 2)


def _q1_symbol_survival() -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for label, path in (("processed", PROCESSED), ("qc", QC)):
        a = ad.read_h5ad(path)
        present = "Igfbp1" in a.var_names
        layer = "counts" if "counts" in a.layers else None
        rows.append(
            {
                "h5ad": label,
                "n_cells": int(a.n_obs),
                "Igfbp1_in_var": present,
                "layer": layer or ".X",
                "pct_expr_all_cells": _pct_expressing(a, "Igfbp1", layer),
            }
        )
    return pd.DataFrame(rows)


def _q2_post_remap() -> pd.DataFrame:
    a = ad.read_h5ad(INTEGRATED)
    mouse = a[a.obs["species"].astype(str) == "mouse"].copy()
    layer = "counts" if "counts" in a.layers else None
    rows = [
        {
            "view": "integrated (all)",
            "n_cells": int(a.n_obs),
            "IGFBP1_in_var": "IGFBP1" in a.var_names,
            "pct_expr": _pct_expressing(a, "IGFBP1", layer),
        },
        {
            "view": "integrated (mouse only)",
            "n_cells": int(mouse.n_obs),
            "IGFBP1_in_var": "IGFBP1" in mouse.var_names,
            "pct_expr": _pct_expressing(mouse, "IGFBP1", layer),
        },
        {
            "view": "integrated (human only)",
            "n_cells": int((a.obs["species"].astype(str) == "human").sum()),
            "IGFBP1_in_var": "IGFBP1" in a.var_names,
            "pct_expr": _pct_expressing(
                a[a.obs["species"].astype(str) == "human"], "IGFBP1", layer
            ),
        },
    ]
    return pd.DataFrame(rows)


def _q3_per_stage_pseudobulk() -> pd.DataFrame:
    a = ad.read_h5ad(QC)
    pb = pseudobulk(a, groupby=["orig.ident", "time"], layer="counts", min_cells=20)
    if "Igfbp1" not in pb.var_names:
        return pd.DataFrame([{"note": "Igfbp1 absent from pseudobulk var set"}])
    j = pb.var_names.get_loc("Igfbp1")
    X = pb.X
    counts = X[:, j].toarray().ravel() if sp.issparse(X) else np.asarray(X[:, j]).ravel()
    lib = X.sum(axis=1)
    lib = np.asarray(lib).ravel()
    cpm = np.where(lib > 0, counts / lib * 1e6, 0.0)
    out = pb.obs.copy()
    out["Igfbp1_counts"] = counts.astype(int)
    out["lib_size"] = lib.astype(int)
    out["Igfbp1_cpm"] = np.round(cpm, 3)
    cols = ["orig.ident", "time", "n_cells", "Igfbp1_counts", "lib_size", "Igfbp1_cpm"]
    return out[cols].sort_values(["time", "orig.ident"]).reset_index(drop=True)


def _verdict(q1: pd.DataFrame, q2: pd.DataFrame, q3: pd.DataFrame) -> str:
    qc_pct = float(q1.loc[q1["h5ad"] == "qc", "pct_expr_all_cells"].iloc[0])
    integ_mouse_pct = float(q2.loc[q2["view"] == "integrated (mouse only)", "pct_expr"].iloc[0])
    if "Igfbp1_counts" in q3.columns:
        stages_with_signal = int((q3["Igfbp1_counts"] > 0).sum())
        total_stages = int(len(q3))
        max_cpm = float(q3["Igfbp1_cpm"].max())
    else:
        stages_with_signal = total_stages = 0
        max_cpm = 0.0

    if qc_pct < 0.5:
        return (
            "**Verdict: real-biology / dataset-capture.** `Igfbp1` is "
            f"expressed in only {qc_pct:.2f} % of mouse cells already at "
            "the QC stage (pre-orthology, mouse symbol). The 0 % in the "
            "integrated joint object reflects the upstream signal, not a "
            "remap bug or HVG dropout. **Action:** IGFBP1 cannot carry a "
            "year-one cross-species claim from this dataset; flag as "
            "human-only-carrying in `configs/markers.yaml::protected_core` "
            "or restrict IGFBP1-based claims to the human side."
        )
    if qc_pct >= 0.5 and integ_mouse_pct < 0.5:
        return (
            "**Verdict: remap / integration pipeline bug.** `Igfbp1` is "
            f"expressed in {qc_pct:.2f} % of mouse cells in the QC h5ad "
            f"but `IGFBP1` is {integ_mouse_pct:.2f} % in the mouse subset "
            "of the integrated h5ad. The signal is being lost between "
            "`results/qc/GSE226417.h5ad` and "
            "`results/integrated/stromal_cross_species.h5ad`. **Action:** "
            "follow-up issue against `src/cell_states/integrate.py::"
            "_remap_mouse_genes` (symbol rename / inner-join / "
            "normalization step). Out of session scope today; log only."
        )
    if stages_with_signal < total_stages:
        return (
            "**Verdict: stage-restricted expression.** `Igfbp1` carries "
            f"non-zero pseudobulk signal in {stages_with_signal} / "
            f"{total_stages} (sample × time) groups (max CPM "
            f"{max_cpm:.1f}). GSE226417 is early-pregnancy only, so the "
            "0 % integrated number is a stage-coverage artifact, not a "
            "biological absence. **Action:** keep IGFBP1 in the "
            "protected core with a documented stage caveat; any IGFBP1 "
            "cross-species claim must condition on matched cycle/"
            "pregnancy stage."
        )
    return (
        "**Verdict: signal present at QC, post-remap and per-stage; "
        "investigate why the integrated `% expressing` metric reads 0.**"
    )


def main() -> int:
    q1 = _q1_symbol_survival()
    q2 = _q2_post_remap()
    q3 = _q3_per_stage_pseudobulk()
    verdict = _verdict(q1, q2, q3)

    print("## IGFBP1 mouse audit")
    print()
    print("Generated by `scripts/diagnose_igfbp1_mouse.py`. Re-run after any")
    print("re-ingest or re-QC of GSE226417.")
    print()
    print("### Q1 — symbol survival (raw mouse symbols, pre-orthology)")
    print()
    print(q1.to_markdown(index=False))
    print()
    print("### Q2 — post-remap survival (integrated joint h5ad)")
    print()
    print(q2.to_markdown(index=False))
    print()
    print("### Q3 — per-(sample × time) pseudobulk on QC mouse h5ad")
    print()
    print(q3.to_markdown(index=False))
    print()
    print(verdict)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
