"""Baseline-priming test (Q4.1).

Hypothesis 1 of the Q4 convergent-evolution question (PLAN.md Q4):
spontaneous-deciduator stroma sits at a **lowered activation threshold**.
If true, the *resting* (unstimulated) stromal population of a spontaneous
deciduator (human) should already score higher on the conserved
``decidual_score`` than the resting stromal population of an induced
deciduator (mouse), and the within-species *priming distance* from
resting to fully-decidualized should be **smaller** in the spontaneous
species.

Operational definitions, derived from atlas ``.obs`` (Q4.1 scope report):

- **Resting baseline** per species: ``cell_type == "stromal_fibroblast"``.
  Holds the celltype constant across species and avoids relying on
  cycle-stage labels that differ between human (cycle phases) and mouse
  (GSE226417 is uniformly early-pregnancy).
- **Decidualized end-state** per species: ``cell_type == "decidual_stromal"``.

Outputs a long table with one row per (species, score) summarising
resting and decidualized means, the within-species Cohen's d
(priming distance), and the between-species Welch's t / Cohen's d at
the resting state.

Decision rule (PLAN.md Q4.1):

- If the human within-species priming distance for ``decidual_score`` is
  meaningfully **smaller** than the mouse priming distance (i.e. human
  resting cells start closer to the decidualized end-state),
  hypothesis 1 is *supported*.
- If the two priming distances are comparable, hypothesis 1 is
  *refuted* and Q4.2 / Q4.3 focus on hypotheses 2 (cis-regulatory
  rewiring) and 4 (stromal-niche pre-priming).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import stats

DEFAULT_RESTING_CELLTYPES: tuple[str, ...] = ("stromal_fibroblast",)
DEFAULT_DECIDUALIZED_CELLTYPES: tuple[str, ...] = (
    "pre_decidual_stromal",
    "decidual_stromal",
    "senescent_decidual",
)


@dataclass(frozen=True)
class BaselinePrimingConfig:
    """Configuration for :func:`baseline_priming`.

    ``resting_celltypes`` / ``decidualized_celltypes`` are sequences
    (``cell_type`` values are pooled). The atlas mouse has only ~11
    fully ``decidual_stromal`` cells but ~3,087 ``pre_decidual_stromal``
    cells, so the decidualized end-state defaults to the union of the
    three decidual-lineage labels (matching the convention used in
    ``src/cell_states/regulators.py``).
    """

    species_key: str = "species"
    cell_type_key: str = "cell_type"
    resting_celltypes: tuple[str, ...] = DEFAULT_RESTING_CELLTYPES
    decidualized_celltypes: tuple[str, ...] = DEFAULT_DECIDUALIZED_CELLTYPES
    min_cells_per_group: int = 20


def _cohens_d(a: np.ndarray, b: np.ndarray) -> float:
    """Pooled-SD Cohen's d for two 1-D arrays. Returns NaN if undefined."""
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return float("nan")
    va, vb = np.var(a, ddof=1), np.var(b, ddof=1)
    pooled = np.sqrt(((na - 1) * va + (nb - 1) * vb) / (na + nb - 2))
    if pooled == 0:
        return float("nan")
    return float((np.mean(b) - np.mean(a)) / pooled)


def baseline_priming(
    adata: AnnData,
    score_cols: list[str],
    config: BaselinePrimingConfig | None = None,
) -> pd.DataFrame:
    """Per-(species, score) resting vs decidualized comparison.

    Returns one row per ``(species, score)`` with:

    - ``n_resting`` / ``n_decidualized``
    - ``resting_mean`` / ``resting_std``
    - ``decidualized_mean`` / ``decidualized_std``
    - ``priming_distance`` — within-species Cohen's d from resting to
      decidualized. Smaller = resting already close to end-state.
    - ``welch_t`` / ``welch_p`` — within-species Welch's t-test
      (resting vs decidualized).

    The caller can pivot on ``species`` to compare priming_distance
    across species per module.
    """
    cfg = config or BaselinePrimingConfig()
    missing = [c for c in (cfg.species_key, cfg.cell_type_key) if c not in adata.obs.columns]
    if missing:
        msg = f"adata.obs missing required columns: {missing}"
        raise ValueError(msg)
    missing_scores = [c for c in score_cols if c not in adata.obs.columns]
    if missing_scores:
        msg = f"adata.obs missing score columns: {missing_scores}"
        raise ValueError(msg)

    species_levels = sorted(adata.obs[cfg.species_key].astype(str).unique())
    rows: list[dict[str, object]] = []

    for species in species_levels:
        sp_mask = adata.obs[cfg.species_key].astype(str) == species
        ct = adata.obs[cfg.cell_type_key].astype(str)
        resting_mask = sp_mask & ct.isin(cfg.resting_celltypes)
        decid_mask = sp_mask & ct.isin(cfg.decidualized_celltypes)
        n_rest = int(resting_mask.sum())
        n_decid = int(decid_mask.sum())

        for score in score_cols:
            vals_rest = adata.obs.loc[resting_mask, score].to_numpy(dtype=float)
            vals_decid = adata.obs.loc[decid_mask, score].to_numpy(dtype=float)
            vals_rest = vals_rest[~np.isnan(vals_rest)]
            vals_decid = vals_decid[~np.isnan(vals_decid)]

            if (
                len(vals_rest) < cfg.min_cells_per_group
                or len(vals_decid) < cfg.min_cells_per_group
            ):
                t_stat, p_val = float("nan"), float("nan")
                d = float("nan")
            else:
                t_stat, p_val = stats.ttest_ind(vals_rest, vals_decid, equal_var=False)
                d = _cohens_d(vals_rest, vals_decid)

            rows.append(
                {
                    "species": species,
                    "score": score,
                    "n_resting": n_rest,
                    "n_decidualized": n_decid,
                    "resting_mean": float(np.mean(vals_rest)) if len(vals_rest) else float("nan"),
                    "resting_std": float(np.std(vals_rest, ddof=1))
                    if len(vals_rest) > 1
                    else float("nan"),
                    "decidualized_mean": float(np.mean(vals_decid))
                    if len(vals_decid)
                    else float("nan"),
                    "decidualized_std": float(np.std(vals_decid, ddof=1))
                    if len(vals_decid) > 1
                    else float("nan"),
                    "priming_distance": d,
                    "welch_t": float(t_stat),
                    "welch_p": float(p_val),
                }
            )

    return pd.DataFrame(rows)


def between_species_resting(
    adata: AnnData,
    score_cols: list[str],
    config: BaselinePrimingConfig | None = None,
) -> pd.DataFrame:
    """Direct between-species comparison at the resting baseline.

    For each score, compares the resting-cell distributions of every
    species pair (one row per (score, species_a, species_b)) using
    Welch's t and Cohen's d. ``d > 0`` means species_b has a higher
    mean than species_a — i.e. ``d_human_vs_mouse > 0`` for
    ``decidual_score`` directly supports hypothesis 1.
    """
    cfg = config or BaselinePrimingConfig()
    species_levels = sorted(adata.obs[cfg.species_key].astype(str).unique())
    if len(species_levels) < 2:
        return pd.DataFrame(
            columns=[
                "score",
                "species_a",
                "species_b",
                "n_a",
                "n_b",
                "mean_a",
                "mean_b",
                "welch_t",
                "welch_p",
                "cohens_d_b_minus_a",
            ]
        )

    ct = adata.obs[cfg.cell_type_key].astype(str)
    resting_mask = ct.isin(cfg.resting_celltypes)
    rows: list[dict[str, object]] = []

    for i, sp_a in enumerate(species_levels):
        for sp_b in species_levels[i + 1 :]:
            mask_a = resting_mask & (adata.obs[cfg.species_key].astype(str) == sp_a)
            mask_b = resting_mask & (adata.obs[cfg.species_key].astype(str) == sp_b)
            for score in score_cols:
                a = adata.obs.loc[mask_a, score].to_numpy(dtype=float)
                b = adata.obs.loc[mask_b, score].to_numpy(dtype=float)
                a = a[~np.isnan(a)]
                b = b[~np.isnan(b)]
                if len(a) < cfg.min_cells_per_group or len(b) < cfg.min_cells_per_group:
                    t_stat, p_val, d = float("nan"), float("nan"), float("nan")
                else:
                    t_stat, p_val = stats.ttest_ind(a, b, equal_var=False)
                    d = _cohens_d(a, b)
                rows.append(
                    {
                        "score": score,
                        "species_a": sp_a,
                        "species_b": sp_b,
                        "n_a": int(len(a)),
                        "n_b": int(len(b)),
                        "mean_a": float(np.mean(a)) if len(a) else float("nan"),
                        "mean_b": float(np.mean(b)) if len(b) else float("nan"),
                        "welch_t": float(t_stat),
                        "welch_p": float(p_val),
                        "cohens_d_b_minus_a": d,
                    }
                )
    return pd.DataFrame(rows)
