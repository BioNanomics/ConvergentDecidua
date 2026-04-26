"""scATAC-seq quality control pipeline.

Filters cells by TSS enrichment and fragment count, then applies
TF-IDF normalization and LSI dimensionality reduction.
"""

from __future__ import annotations

import logging

import anndata as ad
import numpy as np
import scanpy as sc

logger = logging.getLogger(__name__)

DEFAULT_PARAMS = {
    "min_tss_enrichment": 2.0,
    "min_fragments": 1000,
    "max_fragments": 100000,
    "n_components": 50,
}


def qc_scatac(
    adata: ad.AnnData,
    params: dict | None = None,
) -> ad.AnnData:
    """Run scATAC-seq QC pipeline.

    Parameters
    ----------
    adata : ad.AnnData
        Peak × cell AnnData (from scATAC).
    params : dict, optional
        Override default QC parameters.

    Returns
    -------
    ad.AnnData
        Filtered and TF-IDF/LSI transformed AnnData.
    """
    p = {**DEFAULT_PARAMS, **(params or {})}
    n_before = adata.n_obs
    logger.info("Starting scATAC QC: %d cells × %d peaks", adata.n_obs, adata.n_vars)

    # Filter by TSS enrichment if available
    if "tss_enrichment" in adata.obs.columns:
        keep = adata.obs["tss_enrichment"] >= p["min_tss_enrichment"]
        adata = adata[keep].copy()
        logger.info("TSS enrichment filter: %d → %d cells", n_before, adata.n_obs)

    # Filter by fragment count if available
    if "n_fragments" in adata.obs.columns:
        keep = (adata.obs["n_fragments"] >= p["min_fragments"]) & (
            adata.obs["n_fragments"] <= p["max_fragments"]
        )
        adata = adata[keep].copy()
        logger.info("Fragment count filter: → %d cells", adata.n_obs)

    # TF-IDF normalization
    adata = _tfidf(adata)

    # LSI (via PCA on TF-IDF matrix)
    sc.pp.pca(adata, n_comps=min(p["n_components"], adata.n_vars - 1, adata.n_obs - 1))
    # Remove first component (typically correlated with read depth)
    if adata.obsm["X_pca"].shape[1] > 1:
        adata.obsm["X_lsi"] = adata.obsm["X_pca"][:, 1:]
    logger.info("LSI: %d components", adata.obsm.get("X_lsi", adata.obsm["X_pca"]).shape[1])

    adata.obs["qc_passed"] = True
    return adata


def _tfidf(adata: ad.AnnData) -> ad.AnnData:
    """Apply TF-IDF normalization to a peak matrix."""
    import scipy.sparse as sp

    X = adata.X
    if sp.issparse(X):
        X = X.toarray()

    # Term frequency: normalize per cell
    tf = X / (X.sum(axis=1, keepdims=True) + 1e-10)

    # Inverse document frequency
    n_cells = X.shape[0]
    idf = np.log1p(n_cells / (1 + (X > 0).sum(axis=0)))

    tfidf = tf * idf
    adata.X = sp.csr_matrix(tfidf)
    return adata
