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
    processed_dir: Path | None = None,
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
    processed_dir
        Optional ``results/processed/`` directory containing per-dataset
        h5ads. Used to label *why* a canonical marker is missing from the
        joint var set (orthology, inner-join, or HVG selection).
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
    upstream = _scan_processed_var_by_species(processed_dir) if processed_dir else {}
    recovery = _marker_recovery(adata, backbone_path, upstream)
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


def _scan_processed_var_by_species(processed_dir: Path) -> dict[str, set[str]]:
    """Return ``{species: set(uppercase var_names)}`` from per-dataset h5ads.

    Reads each ``*.h5ad`` in ``processed_dir`` in backed mode so we only
    pay for the ``var_names`` + the ``species`` column. Used by the
    drop-audit in :func:`_marker_recovery` to distinguish HVG loss from
    upstream orthology / inner-join loss.
    """
    import contextlib

    import anndata as ad

    by_species: dict[str, set[str]] = {}
    if not processed_dir.exists():
        return by_species
    for path in sorted(processed_dir.glob("*.h5ad")):
        try:
            ad_ = ad.read_h5ad(path, backed="r")
        except Exception as exc:  # pragma: no cover - defensive
            logger.warning("drop-audit: could not read %s (%s)", path, exc)
            continue
        try:
            if "species" not in ad_.obs.columns:
                continue
            species_vals = ad_.obs["species"].astype(str).unique().tolist()
            symbols = {str(v).upper() for v in ad_.var_names}
            for sp in species_vals:
                by_species.setdefault(sp, set()).update(symbols)
        finally:
            with contextlib.suppress(Exception):
                ad_.file.close()
    return by_species


def _classify_drop_reason(
    in_joint: bool,
    symbol: str,
    upstream_by_species: dict[str, set[str]],
) -> str:
    """Label *why* a canonical marker is or isn't in the joint var set.

    Returns one of:
    - ``present`` — in integrated var
    - ``lost_hvg`` — in every species' upstream var set, dropped at HVG step
    - ``lost_inner_join`` — in some but not all species' upstream var sets
    - ``lost_orthology`` — in no species' upstream var set
    - ``unknown`` — no upstream data provided
    """
    if in_joint:
        return "present"
    if not upstream_by_species:
        return "unknown"
    sym = symbol.upper()
    present_in = [sp for sp, vs in upstream_by_species.items() if sym in vs]
    n_species = len(upstream_by_species)
    if len(present_in) == n_species:
        return "lost_hvg"
    if len(present_in) == 0:
        return "lost_orthology"
    return "lost_inner_join"


def _marker_recovery(adata, backbone_path: Path | None, upstream_by_species=None):
    """Per-species fraction expressing each canonical marker + drop-audit.

    The integrated h5ad's ``var_names`` are already harmonized to
    uppercase human symbols (mouse rows were symbol-mapped via the
    ortholog backbone before integration), so we do NOT remap here —
    we just look up the human symbol in ``var_names`` and split by
    species. The ``drop_reason`` column distinguishes HVG loss from
    inner-join loss and orthology loss when ``upstream_by_species`` is
    supplied (built by :func:`_scan_processed_var_by_species`).

    The ``backbone_path`` argument is accepted for API symmetry but
    not used; the joint space carries human symbols.
    """
    import numpy as np
    import pandas as pd

    _ = backbone_path  # unused; symbols are already harmonized
    upstream_by_species = upstream_by_species or {}
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
        in_joint = j is not None
        row["in_joint_var"] = in_joint
        row["drop_reason"] = _classify_drop_reason(in_joint, human_symbol, upstream_by_species)
        if not in_joint:
            for species in species_present:
                row[f"{species}_pct_expr"] = None
            rows.append(row)
            continue
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
        "each canonical decidualization marker. The `drop_reason` "
        "column labels *why* a marker is missing from the joint var "
        "set: `present` — in integrated var; `lost_hvg` — present in "
        "every species' upstream processed h5ad but dropped at HVG "
        "selection (the **fixable** failure mode — see pre-Q3 gate "
        "item A in `PLAN.md` for the `protected_core` carveout); "
        "`lost_inner_join` — present in some species' upstream but "
        "not all; `lost_orthology` — absent from every species' "
        "upstream var set; `unknown` — no upstream data provided.",
        "",
    ]
    if len(recovery):
        lines.append(recovery.to_markdown(index=False))
    else:
        lines.append("_No marker recovery table available._")

    lines.append("")
    return "\n".join(lines)
