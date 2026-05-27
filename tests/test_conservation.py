"""Unit tests for ``src/scoring/conservation.py``."""

from __future__ import annotations

import pandas as pd

from src.scoring.conservation import classify_conservation, summarise_modules


def _fdr_row(module, species, group, dev, fdr):
    return {
        "module": module,
        "species": species,
        "group": group,
        "observed_mean": dev,
        "null_mean": 0.0,
        "fdr": fdr,
    }


def test_classify_conserved_up_and_neutral():
    df = pd.DataFrame(
        [
            _fdr_row("M1", "human", "G", dev=0.5, fdr=0.01),
            _fdr_row("M1", "mouse", "G", dev=0.4, fdr=0.02),
            _fdr_row("M2", "human", "G", dev=0.05, fdr=0.5),
            _fdr_row("M2", "mouse", "G", dev=0.04, fdr=0.6),
        ]
    )
    out = classify_conservation(df, fdr_threshold=0.05)
    classes = dict(zip(out["module"], out["class"], strict=False))
    assert classes["M1"] == "conserved-up"
    assert classes["M2"] == "neutral"


def test_classify_divergent_and_biased():
    df = pd.DataFrame(
        [
            _fdr_row("M1", "human", "G", dev=0.5, fdr=0.01),
            _fdr_row("M1", "mouse", "G", dev=-0.5, fdr=0.01),  # divergent
            _fdr_row("M2", "human", "G", dev=0.5, fdr=0.01),
            _fdr_row("M2", "mouse", "G", dev=0.5, fdr=0.5),  # human-biased
            _fdr_row("M3", "human", "G", dev=-0.5, fdr=0.5),
            _fdr_row("M3", "mouse", "G", dev=-0.5, fdr=0.01),  # mouse-biased-down
        ]
    )
    out = classify_conservation(df, fdr_threshold=0.05)
    classes = dict(zip(out["module"], out["class"], strict=False))
    assert classes["M1"] == "divergent"
    assert classes["M2"] == "human-biased-up"
    assert classes["M3"] == "mouse-biased-down"


def test_summarise_modules_picks_conserved_over_biased():
    df = pd.DataFrame(
        [
            _fdr_row("M1", "human", "G1", dev=0.5, fdr=0.01),
            _fdr_row("M1", "mouse", "G1", dev=0.4, fdr=0.02),  # conserved-up in G1
            _fdr_row("M1", "human", "G2", dev=0.9, fdr=0.001),
            _fdr_row("M1", "mouse", "G2", dev=0.1, fdr=0.5),  # human-biased-up in G2
        ]
    )
    out = classify_conservation(df, fdr_threshold=0.05)
    summary = summarise_modules(out)
    assert summary.iloc[0]["summary_class"] == "conserved-up"
    assert summary.iloc[0]["top_group"] == "G1"
