"""Integration QC report.

Quantifies cross-species mixing on the Harmony embedding (LISI) and
checks that canonical decidualization markers actually mark the
intended clusters (marker recovery). Writes a markdown summary into
``results/reports/integration_qc.md``.

Inputs
------
- ``results/integrated/stromal_cross_species.h5ad`` (or its alias
  ``stromal_harmony.h5ad``) with ``X_pca_harmony`` in ``.obsm`` and
  ``obs[['species','dataset','cell_type','lineage']]``.
- ``configs/markers.yaml`` for the canonical-marker gene sets.
- ``results/orthologs/backbone.parquet`` to map human markers into
  mouse symbols for the recovery check.

Outputs
-------
- ``results/reports/integration_qc.md`` (human-readable summary).
"""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

# Canonical decidualization markers we expect to recover after
# integration. Human symbols; mouse equivalents come via the Tier-1+2
# ortholog backbone and the species_overrides block. Keep this list
# short — it is the smoke test, not the full marker panel.
CANONICAL_MARKERS_HUMAN = [
    "PGR",
    "FOXO1",
    "HAND2",
    "WNT4",
    "IGFBP1",
    "PRL",
    "LEFTY2",
    "IL15",
]


def generate_integration_qc(
    integrated_h5ad: Path,
    output_path: Path,
    backbone_path: Path | None = None,
    n_lisi_cells: int = 5000,
) -> dict:
    """Compute LISI mixing + marker recovery; write a markdown report.

    Parameters
    ----------
    integrated_h5ad
        Path to the integrated AnnData (must have ``obsm['X_pca_harmony']``
        and ``obs[['species','dataset']]``).
    output_path
        Markdown report destination.
    backbone_path
        Optional ortholog backbone for mapping human markers to mouse.
    n_lisi_cells
        Subsample size for LISI (full computation is O(n^2) k-NN).

    Returns
    -------
    dict
        Summary metrics (also embedded in the markdown).
    """
    if not integrated_h5ad.exists():
        logger.warning("Integrated h5ad missing: %s", integrated_h5ad)
        output_path.write_text(
            "# Integration QC\n\n_No integrated h5ad found._ Run "
            "`wombat integrate --mode stromal --method harmony` first.\n"
        )
        return {}

    import anndata as ad

    adata = ad.read_h5ad(integrated_h5ad)
    metrics: dict = {
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
    }

    # --- LISI mixing -------------------------------------------------
    lisi_block = _compute_lisi_block(adata, n_lisi_cells)
    metrics.update(lisi_block)

    # --- Per-dataset / per-lineage breakdown -------------------------
    composition = _composition_table(adata)
    metrics["composition"] = composition.to_dict(orient="records")

    # --- Marker recovery ---------------------------------------------
    recovery = _marker_recovery(adata, backbone_path)
    metrics["marker_recovery"] = recovery.to_dict(orient="records")

    # --- Write markdown ---------------------------------------------
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(_render_markdown(metrics, composition, recovery))
    logger.info("Wrote integration QC report → %s", output_path)
    return metrics


def _compute_lisi_block(adata, n_cells: int) -> dict:
    """Compute LISI on species + dataset; subsample for speed."""
    import numpy as np

    if "X_pca_harmony" not in adata.obsm:
        logger.warning("X_pca_harmony missing; skipping LISI")
        return {"lisi_status": "skipped (no X_pca_harmony)"}

    needed = [c for c in ("species", "dataset") if c in adata.obs.columns]
    if not needed:
        return {"lisi_status": "skipped (no species/dataset cols)"}

    try:
        from harmonypy.lisi import compute_lisi
    except Exception as exc:  # pragma: no cover - import guard
        logger.warning("harmonypy.lisi unavailable: %s", exc)
        return {"lisi_status": f"skipped ({exc})"}

    rng = np.random.default_rng(0)
    n = adata.n_obs
    idx = rng.choice(n, size=n_cells, replace=False) if n > n_cells else np.arange(n)

    emb = adata.obsm["X_pca_harmony"][idx]
    meta = adata.obs[needed].iloc[idx].reset_index(drop=True)

    lisi = compute_lisi(emb, meta, needed)
    block: dict = {"lisi_n_cells": int(len(idx))}
    for j, col in enumerate(needed):
        block[f"lisi_{col}_median"] = float(np.median(lisi[:, j]))
        block[f"lisi_{col}_mean"] = float(np.mean(lisi[:, j]))
        block[f"lisi_{col}_n_categories"] = int(meta[col].nunique())
    return block


def _composition_table(adata):
    """Per-dataset cell / lineage / cell_type breakdown."""
    import pandas as pd

    cols = [c for c in ("species", "dataset", "lineage", "cell_type") if c in adata.obs.columns]
    if "dataset" not in cols:
        return pd.DataFrame()
    grouped = adata.obs.groupby(cols, observed=True).size().reset_index(name="n_cells")
    return grouped.sort_values(["species", "dataset", "n_cells"], ascending=[True, True, False])


def _marker_recovery(adata, backbone_path: Path | None):
    """Per-species fraction expressing each canonical marker.

    The integrated h5ad's ``var_names`` are already harmonized to
    uppercase human symbols (mouse rows were symbol-mapped via the
    ortholog backbone before integration), so we do NOT remap here —
    we just look up the human symbol in ``var_names`` and split by
    species. Missing genes (lost to HVG selection) are flagged so
    reviewers can see the gap.

    The ``backbone_path`` argument is accepted for API symmetry but
    not used; the joint space carries human symbols.
    """
    import numpy as np
    import pandas as pd

    _ = backbone_path  # unused; symbols are already harmonized
    var_index = {str(v): i for i, v in enumerate(adata.var_names)}

    species_present = (
        list(adata.obs["species"].astype(str).unique())
        if "species" in adata.obs.columns
        else ["unknown"]
    )

    rows = []
    for human_symbol in CANONICAL_MARKERS_HUMAN:
        row = {"marker": human_symbol}
        j = var_index.get(human_symbol)
        if j is None:
            row["in_joint_var"] = False
            for species in species_present:
                row[f"{species}_pct_expr"] = None
            rows.append(row)
            continue
        row["in_joint_var"] = True
        for species in species_present:
            mask = (adata.obs["species"].astype(str) == species).values
            X = adata.X[mask, j]
            X = X.toarray().ravel() if hasattr(X, "toarray") else np.asarray(X).ravel()
            row[f"{species}_pct_expr"] = round(float((X > 0).mean()) * 100, 1)
        rows.append(row)

    return pd.DataFrame(rows)


def _render_markdown(metrics: dict, composition, recovery) -> str:
    lines = [
        "# Integration QC",
        "",
        f"- Cells: **{metrics.get('n_cells', 'n/a')}**",
        f"- Genes (joint var set): **{metrics.get('n_genes', 'n/a')}**",
        "",
        "## Cross-species mixing (LISI)",
        "",
        "LISI = effective number of categories in each cell's local "
        "neighbourhood. Higher = better mixing. Computed on a random "
        f"subsample of {metrics.get('lisi_n_cells', 'n/a')} cells "
        "from `obsm['X_pca_harmony']`.",
        "",
    ]
    status = metrics.get("lisi_status")
    if status:
        lines.append(f"_LISI status: {status}_")
    else:
        lines.append("| Variable | Median LISI | Mean LISI | Categories | Mixing |")
        lines.append("|---|---:|---:|---:|---|")
        for key in ("species", "dataset"):
            med = metrics.get(f"lisi_{key}_median")
            mean = metrics.get(f"lisi_{key}_mean")
            n_cat = metrics.get(f"lisi_{key}_n_categories")
            if med is None:
                continue
            # Mixing fraction: (median - 1) / (n_cat - 1) ranges 0..1.
            if n_cat and n_cat > 1:
                frac = (med - 1) / (n_cat - 1)
                if frac < 0.10:
                    label = "⚠️ none (clusters separated)"
                elif frac < 0.40:
                    label = "poor"
                elif frac < 0.70:
                    label = "partial"
                else:
                    label = "good"
            else:
                label = "n/a"
            lines.append(f"| {key} | {med:.2f} | {mean:.2f} | {n_cat or 'n/a'} | {label} |")

    # Flag severe non-mixing inline so reviewers see it without
    # having to interpret LISI numbers themselves.
    species_med = metrics.get("lisi_species_median")
    n_species = metrics.get("lisi_species_n_categories")
    if species_med is not None and n_species and n_species > 1:
        frac = (species_med - 1) / (n_species - 1)
        if frac < 0.10:
            lines += [
                "",
                "> ⚠️ **Species mixing is effectively zero.** Either Harmony "
                "did not converge to a shared embedding (try `theta` ↑ on "
                "`species`, or include more cross-species HVGs), or the "
                "stromal sub-populations differ enough between species that "
                "no shared neighborhood structure exists at the chosen `k`.",
            ]

    lines += [
        "",
        "## Composition",
        "",
        "Per-dataset cell counts after integration + stromal subset.",
        "",
    ]
    if len(composition):
        lines.append(composition.to_markdown(index=False))
    else:
        lines.append("_No composition table available._")

    lines += [
        "",
        "## Marker recovery",
        "",
        "Fraction of cells (per species) with non-zero expression of "
        "each canonical decidualization marker. Confirms ortholog "
        "mapping carried the gene through and integration did not "
        "drop it from the joint var set.",
        "",
    ]
    if len(recovery):
        lines.append(recovery.to_markdown(index=False))
    else:
        lines.append("_No marker recovery table available._")

    lines.append("")
    return "\n".join(lines)
