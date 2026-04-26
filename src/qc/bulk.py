"""Bulk RNA-seq quality control pipeline.

Low-count filtering and normalization for bulk RNA-seq count matrices.
"""

from __future__ import annotations

import logging

import anndata as ad
import numpy as np

logger = logging.getLogger(__name__)

DEFAULT_PARAMS = {
    "min_counts_per_sample": 1000,
    "min_counts_per_gene": 10,
    "min_samples_per_gene": 2,
}


def qc_bulk(
    adata: ad.AnnData,
    params: dict | None = None,
) -> ad.AnnData:
    """Run bulk RNA-seq QC pipeline.

    Parameters
    ----------
    adata : ad.AnnData
        Bulk count AnnData (samples × genes).
    params : dict, optional
        Override default QC parameters.

    Returns
    -------
    ad.AnnData
        Filtered and normalized AnnData.
    """
    p = {**DEFAULT_PARAMS, **(params or {})}
    logger.info("Starting bulk RNA QC: %d samples × %d genes", adata.n_obs, adata.n_vars)

    # Ensure X is dense for bulk (small matrix) or at least csr
    import scipy.sparse

    if scipy.sparse.issparse(adata.X):
        adata.X = np.asarray(adata.X.todense())

    # Filter low-count samples
    sample_counts = np.asarray(adata.X.sum(axis=1)).flatten()
    keep_samples = sample_counts >= p["min_counts_per_sample"]
    adata = adata[keep_samples].copy()
    logger.info("Sample filter: %d samples retained", adata.n_obs)

    # Filter low-count genes
    gene_counts = np.asarray(adata.X.sum(axis=0)).flatten()
    gene_detected = np.asarray((adata.X > 0).sum(axis=0)).flatten()
    keep_genes = (gene_counts >= p["min_counts_per_gene"]) & (
        gene_detected >= p["min_samples_per_gene"]
    )
    adata = adata[:, keep_genes].copy()
    logger.info("Gene filter: %d genes retained", adata.n_vars)

    # Store raw counts
    adata.layers["counts"] = adata.X.copy()

    # CPM normalization + log
    total = np.asarray(adata.X.sum(axis=1)).flatten()
    adata.X = np.log1p(adata.X / total[:, None] * 1e6)

    adata.obs["qc_passed"] = True
    return adata
