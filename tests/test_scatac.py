"""Unit tests for src/qc/scatac.py.

These tests run on tiny synthetic matrices — no real scATAC data needed.
They lock in two Q2.5 fixes:

1. ``_tfidf`` stays sparse end-to-end (real scATAC matrices are
   ~10^4 cells × 10^5 peaks; a stray ``.toarray()`` OOMs the box).
2. ``gene_activity`` correctly aggregates peak counts into a per-gene
   window around each TSS, respecting strand.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from src.qc.scatac import _tfidf, gene_activity


def _make_peak_adata(n_cells: int = 6, n_peaks: int = 5):
    rng = np.random.default_rng(0)
    X = sp.csr_matrix(rng.integers(0, 4, size=(n_cells, n_peaks)))
    var = pd.DataFrame(
        {
            "chrom": ["chr1"] * n_peaks,
            "start": [100, 500, 1500, 3000, 5000],
            "end":   [200, 600, 1600, 3100, 5100],
        },
        index=[f"peak{i}" for i in range(n_peaks)],
    )
    return ad.AnnData(X=X, var=var)


def test_tfidf_stays_sparse():
    a = _make_peak_adata()
    out = _tfidf(a)
    assert sp.issparse(out.X), "TF-IDF output must stay sparse"
    assert out.X.shape == (6, 5)
    # Per-cell row sums of TF should each be ≤ idf_max (sanity: finite,
    # non-negative, no NaN).
    arr = out.X.toarray()
    assert np.all(np.isfinite(arr))
    assert np.all(arr >= 0)


def test_tfidf_handles_zero_row():
    a = _make_peak_adata()
    # Force one all-zero cell — the legacy implementation divided by
    # row_sum + 1e-10 and produced near-zero floats; the new one uses
    # row_sums[row_sums == 0] = 1.0 and should produce an exact zero row.
    X = a.X.toarray()
    X[2, :] = 0
    a.X = sp.csr_matrix(X)
    out = _tfidf(a)
    assert sp.issparse(out.X)
    row2 = out.X.toarray()[2, :]
    assert np.all(row2 == 0)


def test_gene_activity_window_and_strand():
    peaks = _make_peak_adata(n_cells=4, n_peaks=5)
    # Peak midpoints: 150, 550, 1550, 3050, 5050.
    genes = pd.DataFrame(
        {
            "gene_symbol": ["GENE_A", "GENE_B"],
            "chrom": ["chr1", "chr1"],
            # GENE_A TSS=1000 on + strand, window [-2000, 0] => [-1000, 1000]
            #   → peaks at 150 and 550 fall in window.
            # GENE_B TSS=4000 on - strand, window flips so it's
            #   [TSS - downstream, TSS + upstream] => [4000, 6000]
            #   → peaks at 5050 fall in window.
            "tss": [1000, 4000],
            "strand": ["+", "-"],
        }
    )
    out = gene_activity(peaks, genes, upstream=2000, downstream=0)
    assert out.n_obs == 4
    assert list(out.var.index) == ["GENE_A", "GENE_B"]

    # Expected GENE_A column = sum of peak0 (mid 150) + peak1 (mid 550)
    peak_mat = peaks.X.toarray()
    expected_a = peak_mat[:, 0] + peak_mat[:, 1]
    # Expected GENE_B column = peak4 (mid 5050)
    expected_b = peak_mat[:, 4]

    got = out.X.toarray()
    np.testing.assert_array_equal(got[:, 0], expected_a)
    np.testing.assert_array_equal(got[:, 1], expected_b)


def test_gene_activity_missing_columns_raises():
    peaks = _make_peak_adata()
    bad_genes = pd.DataFrame({"gene_symbol": ["X"], "chrom": ["chr1"]})
    import pytest

    with pytest.raises(ValueError, match="gene_coords missing"):
        gene_activity(peaks, bad_genes)
