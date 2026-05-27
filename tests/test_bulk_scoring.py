"""Unit tests for ``src/scoring/bulk.py`` time-axis parser and monotonicity."""

from __future__ import annotations

import anndata as ad
import numpy as np
import pandas as pd
import pytest

from src.scoring.bulk import monotonicity, parse_time_axis


def test_parse_time_axis_known_labels():
    s = parse_time_axis(["Control", "Day1", "Day2", "Day5", "D0", "day_3"])
    assert s.tolist() == [0.0, 1.0, 2.0, 5.0, 0.0, 3.0]


def test_parse_time_axis_unknown_label_is_nan():
    s = parse_time_axis(["foo", "Day1"])
    assert np.isnan(s.iloc[0])
    assert s.iloc[1] == 1.0


def _toy_scored_adata(scores: dict[str, list[float]]) -> ad.AnnData:
    n = len(next(iter(scores.values())))
    a = ad.AnnData(X=np.zeros((n, 1)))
    a.obs_names = [f"Day{i}" for i in range(n)]
    a.obs = pd.DataFrame(scores, index=a.obs_names)
    a.obs["time"] = parse_time_axis(list(a.obs_names)).values
    return a


def test_monotonicity_perfect_increase():
    a = _toy_scored_adata({"decidual_score": [0.1, 0.5, 0.9, 1.3, 1.7, 2.1]})
    table = monotonicity(a, ["decidual_score"])
    row = table.loc["decidual_score"]
    assert row["rho"] == pytest.approx(1.0)
    assert bool(row["monotonic"]) is True
    assert row["n_samples"] == 6


def test_monotonicity_below_threshold():
    a = _toy_scored_adata({"decidual_score": [0.1, 0.2, 0.1, 0.3, 0.2, 0.4]})
    row = monotonicity(a, ["decidual_score"]).loc["decidual_score"]
    # Modest positive trend but not at the |rho|>=0.7 + p<0.05 floor.
    assert bool(row["monotonic"]) is False


def test_monotonicity_requires_time():
    a = ad.AnnData(X=np.zeros((3, 1)))
    a.obs = pd.DataFrame({"decidual_score": [0.0, 0.5, 1.0]}, index=["a", "b", "c"])
    with pytest.raises(ValueError, match="time"):
        monotonicity(a, ["decidual_score"])
