"""scRNA-seq quality control pipeline.

Filters cells and genes, detects doublets, normalizes, and selects
highly variable genes. QC metrics are stored in .obs.
"""

from __future__ import annotations

import logging

import anndata as ad
import scanpy as sc

logger = logging.getLogger(__name__)

# Default QC parameters (per AI Phase Plan)
DEFAULT_PARAMS = {
    "human": {
        "min_genes": 500,
        "max_genes": 8000,
        "max_pct_mito": 15.0,
        "min_cells": 3,
    },
    "mouse": {
        "min_genes": 500,
        "max_genes": 7000,
        "max_pct_mito": 15.0,
        "min_cells": 3,
    },
}


def qc_scrna(
    adata: ad.AnnData,
    species: str = "human",
    params: dict | None = None,
    *,
    detect_doublets: bool = True,
) -> ad.AnnData:
    """Run scRNA-seq QC pipeline.

    Parameters
    ----------
    adata : ad.AnnData
        Raw count AnnData.
    species : str
        Species name for parameter defaults.
    params : dict, optional
        Override default QC parameters.
    detect_doublets : bool
        Whether to run doublet detection.

    Returns
    -------
    ad.AnnData
        Filtered, normalized AnnData with QC metrics in .obs.
    """
    p = {**DEFAULT_PARAMS.get(species, DEFAULT_PARAMS["human"]), **(params or {})}
    n_before = adata.n_obs
    logger.info("Starting scRNA QC: %d cells × %d genes", adata.n_obs, adata.n_vars)

    # Compute QC metrics
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    # Filter cells
    keep = (
        (adata.obs["n_genes_by_counts"] >= p["min_genes"])
        & (adata.obs["n_genes_by_counts"] <= p["max_genes"])
        & (adata.obs["pct_counts_mt"] <= p["max_pct_mito"])
    )
    adata = adata[keep].copy()
    logger.info("Cell filter: %d → %d cells", n_before, adata.n_obs)

    # Filter genes
    sc.pp.filter_genes(adata, min_cells=p["min_cells"])
    logger.info("Gene filter: %d genes retained", adata.n_vars)

    # Doublet detection
    if detect_doublets:
        adata = _detect_doublets(adata)

    # Store raw counts before normalization
    adata.layers["counts"] = adata.X.copy()

    # Normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # HVG selection — use seurat_v3 if skmisc available, fall back to cell_ranger
    try:
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="seurat_v3", layer="counts")
    except (ImportError, ValueError):
        logger.warning("seurat_v3 HVG failed, falling back to cell_ranger flavor")
        sc.pp.highly_variable_genes(adata, n_top_genes=3000, flavor="cell_ranger")
    logger.info("Selected %d highly variable genes", adata.var["highly_variable"].sum())

    adata.obs["qc_passed"] = True
    return adata


def _detect_doublets(adata: ad.AnnData) -> ad.AnnData:
    """Run doublet detection using scrublet-like approach via scanpy."""
    try:
        sc.pp.scrublet(adata, verbose=False)
        n_doublets = adata.obs.get("predicted_doublet", pd.Series()).sum()
        logger.info("Doublet detection: %d predicted doublets", n_doublets)

        if "predicted_doublet" in adata.obs.columns:
            adata = adata[~adata.obs["predicted_doublet"]].copy()
            logger.info("Removed doublets: %d cells remaining", adata.n_obs)
    except Exception as exc:
        logger.warning("Doublet detection failed (continuing without): %s", exc)
        adata.obs["doublet_score"] = 0.0
        adata.obs["predicted_doublet"] = False

    return adata


# Needed for the fallback in _detect_doublets
import pandas as pd  # noqa: E402
