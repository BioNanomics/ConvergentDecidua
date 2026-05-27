"""Unit tests for ``src/scoring/null.py``."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest
import scipy.sparse

from src.scoring.engine import score_all_modules
from src.scoring.null import NullConfig, _bh, score_with_null


def _toy_adata(seed: int = 0) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    n_cells = 80
    n_genes = 200
    # Build log-norm expression: most genes ~ N(0,1), a "decidual" set in
    # cell_state=B inflated by +3, "noise" set never elevated.
    X = rng.normal(0.0, 1.0, size=(n_cells, n_genes)).astype(np.float32)
    decidual_idx = np.arange(0, 10)
    cell_state = np.array(["A"] * 40 + ["B"] * 40)
    # Inflate decidual genes in cell_state B (loop because fancy-index slices
    # on csr-bound arrays don't write back through chained indexing).
    for i in np.where(cell_state == "B")[0]:
        X[i, decidual_idx] += 3.0
    var_names = [f"g{i}" for i in range(n_genes)]
    obs = pd.DataFrame(
        {
            "species": "human",
            "cell_state": cell_state,
        },
        index=[f"c{i}" for i in range(n_cells)],
    )
    a = ad.AnnData(X=scipy.sparse.csr_matrix(X), obs=obs)
    a.var_names = var_names
    return a


def test_bh_known_input():
    # BH: q_i = min_{k>=i} ( p_k * n / k ), clipped monotone.
    p = np.array([0.001, 0.008, 0.039, 0.041, 0.042])
    q = _bh(p)
    expected = np.array([0.005, 0.020, 0.042, 0.042, 0.042])
    np.testing.assert_allclose(q, expected, atol=1e-6)


def test_bh_monotone_and_bounded():
    rng = np.random.default_rng(0)
    p = rng.uniform(size=50)
    q = _bh(p)
    assert q.min() >= 0.0
    assert q.max() <= 1.0
    # values in p-rank order must be non-decreasing
    order = np.argsort(p)
    assert np.all(np.diff(q[order]) >= -1e-12)


def test_score_with_null_signal_vs_noise():
    adata = _toy_adata()
    gene_sets = {
        "decidual_score": [f"g{i}" for i in range(10)],
        "noise_score": [f"g{i}" for i in range(20, 30)],
    }
    adata = score_all_modules(adata, gene_sets, species="human", backbone_path=None)
    out = score_with_null(
        adata,
        gene_sets,
        species_to_backbone={"human": None},
        config=NullConfig(n_permutations=100, seed=7),
    )
    # Real signal in cell_state B for decidual_score should be highly significant;
    # noise score should not be.
    dec_b = out[(out["module"] == "decidual_score") & (out["group"] == "B")].iloc[0]
    noise_a = out[(out["module"] == "noise_score") & (out["group"] == "A")].iloc[0]
    assert dec_b["fdr"] < 0.05
    assert noise_a["fdr"] > 0.1


def test_score_with_null_requires_scored_obs():
    adata = _toy_adata()
    with pytest.raises(ValueError, match="missing"):
        score_with_null(
            adata,
            {"decidual_score": [f"g{i}" for i in range(10)]},
            species_to_backbone={"human": None},
            config=NullConfig(n_permutations=10),
        )
