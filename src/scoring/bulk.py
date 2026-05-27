"""Score bulk RNA-seq datasets with the same modules used on scRNA.

Bulk samples are tiny (n_samples in the single digits), so the same
``scanpy.tl.score_genes`` engine works fine — we just need a stable
sample-level time axis to compute monotonicity of e.g.
``decidual_score`` against the experimental day.

Q3.1 prerequisite: validates that the scoring modules are not single-
cell-specific. Monotonic ``decidual_score`` vs. day in an in-vitro
decidualization time course is a sanity floor.
"""

from __future__ import annotations

import logging
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

from src.scoring.engine import score_all_modules
from src.scoring.gene_sets import load_score_gene_sets

logger = logging.getLogger(__name__)


def parse_time_axis(obs_names: list[str]) -> pd.Series:
    """Extract a numeric time axis from bulk sample labels.

    Recognises the pattern ``Control``/``D0``/``Day0`` → 0 and
    ``Day<n>`` / ``D<n>`` → ``n``. Returns NaN for unrecognised labels;
    the caller decides whether that is fatal.
    """
    values: list[float] = []
    for name in obs_names:
        s = str(name).strip()
        if s.lower() in ("control", "ctrl", "vehicle", "d0", "day0"):
            values.append(0.0)
            continue
        m = re.search(r"d(?:ay)?[_\-\s]?(\d+(?:\.\d+)?)", s, re.IGNORECASE)
        if m:
            values.append(float(m.group(1)))
            continue
        values.append(float("nan"))
    return pd.Series(values, index=list(obs_names), name="time", dtype=float)


def score_bulk(
    adata: ad.AnnData,
    species: str,
    backbone_path: Path | None,
    gene_sets: dict[str, list[str]] | None = None,
) -> ad.AnnData:
    """Apply ``score_all_modules`` to a bulk AnnData (samples × genes).

    Adds a numeric ``time`` column to ``.obs`` derived from
    ``obs_names`` via :func:`parse_time_axis`.
    """
    if gene_sets is None:
        gene_sets = load_score_gene_sets()
    adata = score_all_modules(adata, gene_sets, species=species, backbone_path=backbone_path)
    adata.obs["time"] = parse_time_axis(list(adata.obs_names)).values
    return adata


def monotonicity(adata: ad.AnnData, score_cols: list[str]) -> pd.DataFrame:
    """Spearman rank correlation between ``time`` and each score column.

    Returns a DataFrame indexed by score with columns
    ``rho``, ``pval``, ``n_samples``, ``monotonic``. ``monotonic`` is
    True when |rho| ≥ 0.7 and pval < 0.05 (heuristic floor; small-n
    bulk has limited power so the rho floor is the primary signal).
    """
    if "time" not in adata.obs.columns:
        msg = "adata.obs['time'] missing — call score_bulk() first"
        raise ValueError(msg)

    time = adata.obs["time"].astype(float).values
    valid = ~np.isnan(time)
    if valid.sum() < 3:
        msg = f"Need ≥3 samples with a valid time axis, got {int(valid.sum())}"
        raise ValueError(msg)

    rows: list[dict[str, float | bool | str]] = []
    for col in score_cols:
        if col not in adata.obs.columns:
            continue
        y = adata.obs[col].astype(float).values
        ok = valid & ~np.isnan(y)
        if ok.sum() < 3:
            rows.append(
                {
                    "score": col,
                    "rho": float("nan"),
                    "pval": float("nan"),
                    "n_samples": int(ok.sum()),
                    "monotonic": False,
                }
            )
            continue
        rho, pval = stats.spearmanr(time[ok], y[ok])
        rows.append(
            {
                "score": col,
                "rho": round(float(rho), 4),
                "pval": round(float(pval), 4),
                "n_samples": int(ok.sum()),
                "monotonic": bool(abs(rho) >= 0.7 and pval < 0.05),
            }
        )
    return pd.DataFrame(rows).set_index("score")
