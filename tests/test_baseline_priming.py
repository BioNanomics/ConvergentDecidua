"""Unit tests for src/scoring/baseline_priming.py (Q4.1)."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from anndata import AnnData

from src.scoring.baseline_priming import (
    BaselinePrimingConfig,
    baseline_priming,
    between_species_resting,
)


def _make_atlas(rng: np.random.Generator) -> AnnData:
    """Synthetic atlas: human stromal_fibroblast starts close to the
    decidualized end-state (small priming distance, supports H1);
    mouse stromal_fibroblast starts far from it (large priming
    distance, does not support H1)."""
    n_human_rest, n_human_decid = 300, 100
    n_mouse_rest, n_mouse_decid = 400, 80

    decidual_score = np.concatenate(
        [
            rng.normal(loc=0.6, scale=0.2, size=n_human_rest),  # human resting: primed
            rng.normal(loc=1.0, scale=0.2, size=n_human_decid),  # human decid
            rng.normal(loc=0.0, scale=0.2, size=n_mouse_rest),  # mouse resting: flat
            rng.normal(loc=1.0, scale=0.2, size=n_mouse_decid),  # mouse decid
        ]
    )
    # An unrelated score with the SAME priming distance in both species.
    other_score = np.concatenate(
        [
            rng.normal(loc=0.0, scale=0.2, size=n_human_rest),
            rng.normal(loc=1.0, scale=0.2, size=n_human_decid),
            rng.normal(loc=0.0, scale=0.2, size=n_mouse_rest),
            rng.normal(loc=1.0, scale=0.2, size=n_mouse_decid),
        ]
    )
    obs = pd.DataFrame(
        {
            "species": ["human"] * (n_human_rest + n_human_decid)
            + ["mouse"] * (n_mouse_rest + n_mouse_decid),
            "cell_type": ["stromal_fibroblast"] * n_human_rest
            + ["decidual_stromal"] * n_human_decid
            + ["stromal_fibroblast"] * n_mouse_rest
            + ["decidual_stromal"] * n_mouse_decid,
            "decidual_score": decidual_score,
            "other_score": other_score,
        }
    )
    n = len(obs)
    return AnnData(X=np.zeros((n, 1)), obs=obs)


def test_baseline_priming_recovers_planted_signal():
    adata = _make_atlas(np.random.default_rng(0))
    out = baseline_priming(adata, score_cols=["decidual_score", "other_score"])
    pivot = out.set_index(["species", "score"])

    human_d = pivot.loc[("human", "decidual_score"), "priming_distance"]
    mouse_d = pivot.loc[("mouse", "decidual_score"), "priming_distance"]
    # decidual_score: human priming distance should be markedly smaller
    # than mouse (~2 sd vs ~5 sd in the synthetic data).
    assert human_d < mouse_d
    assert mouse_d - human_d > 1.0

    # other_score: both species have ~5 sd priming distance.
    human_o = pivot.loc[("human", "other_score"), "priming_distance"]
    mouse_o = pivot.loc[("mouse", "other_score"), "priming_distance"]
    assert abs(human_o - mouse_o) < 0.5

    # Means + counts.
    assert pivot.loc[("human", "decidual_score"), "n_resting"] == 300
    assert pivot.loc[("mouse", "decidual_score"), "n_resting"] == 400
    assert pivot.loc[("human", "decidual_score"), "resting_mean"] > 0.4


def test_baseline_priming_min_cells_floor():
    adata = _make_atlas(np.random.default_rng(0))
    out = baseline_priming(
        adata,
        score_cols=["decidual_score"],
        config=BaselinePrimingConfig(min_cells_per_group=10_000),
    )
    # No species clears the floor → all NaN.
    assert out["priming_distance"].isna().all()
    assert out["welch_p"].isna().all()


def test_between_species_resting_recovers_human_higher():
    adata = _make_atlas(np.random.default_rng(0))
    out = between_species_resting(adata, score_cols=["decidual_score"])
    assert len(out) == 1
    row = out.iloc[0]
    # alphabetical sort → species_a = "human", species_b = "mouse".
    assert (row["species_a"], row["species_b"]) == ("human", "mouse")
    # mouse - human ≈ 0.0 - 0.6 = -0.6  → cohens_d_b_minus_a is negative
    # (mouse resting LOWER than human resting → supports H1).
    assert row["cohens_d_b_minus_a"] < -1.0
    assert row["welch_p"] < 1e-10


def test_baseline_priming_rejects_missing_columns():
    adata = AnnData(X=np.zeros((3, 1)), obs=pd.DataFrame({"species": ["human"] * 3}))
    with pytest.raises(ValueError, match="missing required columns"):
        baseline_priming(adata, score_cols=["decidual_score"])

    obs = pd.DataFrame({"species": ["human"] * 3, "cell_type": ["stromal_fibroblast"] * 3})
    adata2 = AnnData(X=np.zeros((3, 1)), obs=obs)
    with pytest.raises(ValueError, match="missing score columns"):
        baseline_priming(adata2, score_cols=["decidual_score"])
