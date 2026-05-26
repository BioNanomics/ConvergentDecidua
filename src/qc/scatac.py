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
    """Apply TF-IDF normalization to a peak matrix (sparse-safe).

    Real scATAC matrices are O(10^4 cells × 10^5 peaks); a dense
    toarray() blows up to tens of GB. This implementation stays in
    CSR throughout and uses scipy.sparse element-wise multiplication
    with row/column scaling vectors.
    """
    import scipy.sparse as sp

    X = adata.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)

    # Term frequency: divide each row by its sum.
    row_sums = np.asarray(X.sum(axis=1)).ravel()
    row_sums[row_sums == 0] = 1.0
    tf = sp.diags(1.0 / row_sums) @ X  # still sparse

    # Inverse document frequency on the binary (peak-open) signal.
    n_cells = X.shape[0]
    # (X > 0).sum(axis=0) — sparse-safe via .getnnz on the boolean view
    peak_open = np.asarray((X != 0).sum(axis=0)).ravel().astype(float)
    idf = np.log1p(n_cells / (1.0 + peak_open))

    # Scale columns by idf.
    tfidf = tf @ sp.diags(idf)
    adata.X = tfidf.tocsr()
    return adata


def gene_activity(
    adata_peaks: ad.AnnData,
    gene_coords,
    upstream: int = 2000,
    downstream: int = 0,
) -> ad.AnnData:
    """Aggregate peak signal into a cells × genes activity matrix.

    Signac-style: for each gene, sum the counts of all peaks whose
    midpoint falls within [TSS - ``upstream``, TSS + ``downstream``]
    (or, for genes on the minus strand, the analogous window on the
    other side of the TSS).

    Parameters
    ----------
    adata_peaks : ad.AnnData
        Cells × peaks AnnData. ``var`` must carry ``chrom``, ``start``,
        ``end`` columns (peak coordinates in 0-based half-open format).
    gene_coords : pandas.DataFrame
        Gene coordinates with columns ``gene_symbol``, ``chrom``,
        ``tss``, ``strand``. One row per gene.
    upstream, downstream : int
        Window around each TSS to aggregate peaks over.

    Returns
    -------
    ad.AnnData
        Cells × genes AnnData with summed peak counts as ``.X``. Same
        obs as the input; var indexed by ``gene_symbol``.
    """
    import pandas as pd
    import scipy.sparse as sp

    required_var = {"chrom", "start", "end"}
    missing = required_var - set(adata_peaks.var.columns)
    if missing:
        raise ValueError(f"adata_peaks.var missing required columns: {sorted(missing)}")

    required_gene = {"gene_symbol", "chrom", "tss", "strand"}
    missing = required_gene - set(gene_coords.columns)
    if missing:
        raise ValueError(f"gene_coords missing required columns: {sorted(missing)}")

    peak_var = adata_peaks.var.copy()
    # Track positional index so we can slice X.tocsc()[:, idx] later.
    peak_var["__pos__"] = np.arange(len(peak_var))
    peak_var["mid"] = (peak_var["start"].astype(int) + peak_var["end"].astype(int)) // 2

    # Build per-chrom indexes for fast lookup
    peaks_by_chrom: dict[str, pd.DataFrame] = {
        chrom: grp.sort_values("mid").reset_index(drop=True)
        for chrom, grp in peak_var.groupby("chrom")
    }

    n_cells = adata_peaks.n_obs
    n_genes = len(gene_coords)

    X = adata_peaks.X
    if not sp.issparse(X):
        X = sp.csr_matrix(X)
    X = X.tocsc()  # column slicing for peak indices

    rows, cols, data = [], [], []
    for g_idx, gene in enumerate(gene_coords.itertuples(index=False)):
        chrom = str(gene.chrom)
        peaks = peaks_by_chrom.get(chrom)
        if peaks is None or len(peaks) == 0:
            continue
        tss = int(gene.tss)
        if str(gene.strand) == "-":
            lo, hi = tss - downstream, tss + upstream
        else:
            lo, hi = tss - upstream, tss + downstream
        if lo > hi:
            lo, hi = hi, lo
        hit_mask = (peaks["mid"] >= lo) & (peaks["mid"] <= hi)
        if not hit_mask.any():
            continue
        peak_idx = peaks.loc[hit_mask, "__pos__"].astype(int).values
        sub = X[:, peak_idx]
        col_sum = np.asarray(sub.sum(axis=1)).ravel()
        nz = col_sum.nonzero()[0]
        if len(nz) == 0:
            continue
        rows.extend(nz.tolist())
        cols.extend([g_idx] * len(nz))
        data.extend(col_sum[nz].tolist())

    activity = sp.coo_matrix((data, (rows, cols)), shape=(n_cells, n_genes)).tocsr()
    var = pd.DataFrame(index=gene_coords["gene_symbol"].astype(str).values)
    out = ad.AnnData(X=activity, obs=adata_peaks.obs.copy(), var=var)
    logger.info(
        "Gene activity: %d cells × %d genes, %d non-zero entries",
        out.n_obs,
        out.n_vars,
        activity.nnz,
    )
    return out
