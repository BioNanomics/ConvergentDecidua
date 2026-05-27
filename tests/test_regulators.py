"""Unit tests for ``src/cell_states/regulators.py``."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse

from src.cell_states.regulators import (
    RegulatorConfig,
    rank_regulators,
    split_regulators,
)


def _toy_lineage_adata(seed: int = 0) -> ad.AnnData:
    rng = np.random.default_rng(seed)
    n_cells = 200
    n_genes = 30
    X = rng.normal(0.0, 0.2, size=(n_cells, n_genes)).astype(np.float32).clip(0)
    score = rng.normal(0.0, 1.0, size=n_cells).astype(np.float32)
    # Plant a strong-correlation TF and a noise TF.
    X[:, 0] = score * 1.5 + rng.normal(0.0, 0.1, size=n_cells)  # TF_HIGH
    X[:, 1] = rng.normal(0.0, 1.0, size=n_cells)  # TF_NOISE
    var_names = ["TF_HIGH", "TF_NOISE"] + [f"g{i}" for i in range(2, n_genes)]
    obs = pd.DataFrame(
        {
            "species": np.array(["human"] * 100 + ["mouse"] * 100),
            "cell_type": np.array(["decidual_stromal"] * n_cells),
            "decidual_score": score,
        },
        index=[f"c{i}" for i in range(n_cells)],
    )
    a = ad.AnnData(X=scipy.sparse.csr_matrix(X), obs=obs)
    a.var_names = var_names
    return a


def test_rank_regulators_finds_signal():
    adata = _toy_lineage_adata()
    tfs = ["TF_HIGH", "TF_NOISE"]
    ranked = rank_regulators(adata, tfs, config=RegulatorConfig(min_nonzero_fraction=0.0))
    by_tf = ranked.set_index("tf")
    # TF_HIGH must rank above TF_NOISE in both species and have higher |rho|.
    assert by_tf.loc["TF_HIGH", "mean_rank"] < by_tf.loc["TF_NOISE", "mean_rank"]
    assert abs(by_tf.loc["TF_HIGH", "human_rho"]) > abs(by_tf.loc["TF_NOISE", "human_rho"])


def test_split_regulators_caps_correctly():
    df = pd.DataFrame(
        {
            "tf": [f"T{i}" for i in range(30)],
            "human_rho": np.linspace(0.9, 0.1, 30),
            "mouse_rho": np.linspace(0.9, 0.1, 30),
            "human_rank": np.arange(1, 31),
            "mouse_rank": np.arange(1, 31),
            "mean_rank": np.arange(1, 31),
            "rank_gap": np.zeros(30),
        }
    )
    splits = split_regulators(df, cap=10)
    assert len(splits["conserved"]) == 10
    assert splits["conserved"].iloc[0]["tf"] == "T0"
