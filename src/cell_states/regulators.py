"""Candidate regulator shortlists for the decidualization axis (Q3.4).

Reads the integrated/scored atlas, the curated Lambert-2018 human TF
list, and the per-cell ``decidual_score`` and ranks TFs by how strongly
their expression correlates with the decidualization axis in each
species.

Two lists are emitted:

- **Conserved regulators** — TFs ranked high in BOTH species (low mean
  of the two ranks). Capped at ~25.
- **Divergent regulators** — TFs with the largest rank gap between
  species (separately for human-only and mouse-only directions). Capped
  at ~25 each direction.

Scope: cells with ``cell_type in {"pre_decidual_stromal",
"decidual_stromal", "senescent_decidual"}`` (the decidual lineage).
Correlation is Spearman of cell-level TF expression vs.
``decidual_score`` within that lineage in each species.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from scipy import stats

logger = logging.getLogger(__name__)

DEFAULT_DECIDUAL_LINEAGE = (
    "pre_decidual_stromal",
    "decidual_stromal",
    "senescent_decidual",
)

DEFAULT_TF_LIST = (
    Path(__file__).resolve().parents[2] / "configs/reference/lambert2018_human_TFs.txt"
)


@dataclass(frozen=True)
class RegulatorConfig:
    score_col: str = "decidual_score"
    species_key: str = "species"
    cell_type_key: str = "cell_type"
    min_cells_per_species: int = 30
    min_nonzero_fraction: float = 0.05  # drop TFs detected in <5% of cells
    cap: int = 25


def load_tf_list(path: Path | str | None = None) -> list[str]:
    """Load a one-symbol-per-line TF list (Lambert 2018 by default)."""
    p = Path(path) if path else DEFAULT_TF_LIST
    if not p.exists():
        msg = f"TF list not found at {p}"
        raise FileNotFoundError(msg)
    return [ln.strip() for ln in p.read_text().splitlines() if ln.strip()]


def rank_regulators(
    adata: ad.AnnData,
    tf_list: list[str],
    config: RegulatorConfig | None = None,
    lineage: tuple[str, ...] = DEFAULT_DECIDUAL_LINEAGE,
) -> pd.DataFrame:
    """Per-species Spearman correlation between TF expression and decidual_score.

    Returns a wide table indexed by TF symbol with columns
    ``{species}_rho``, ``{species}_pval``, ``{species}_rank``,
    plus ``mean_rank`` and ``rank_gap``.
    """
    cfg = config or RegulatorConfig()
    if cfg.score_col not in adata.obs.columns:
        msg = f"adata.obs['{cfg.score_col}'] missing"
        raise ValueError(msg)

    in_lineage = adata.obs[cfg.cell_type_key].astype(str).isin(lineage).values
    sub = adata[in_lineage].copy()
    if sub.n_obs == 0:
        msg = f"No cells in lineage {lineage}"
        raise ValueError(msg)
    logger.info("Decidual-lineage cells: %d / %d", sub.n_obs, adata.n_obs)

    species_levels = sorted(map(str, sub.obs[cfg.species_key].unique()))
    tfs_present = [tf for tf in tf_list if tf in sub.var_names]
    logger.info("TFs present in atlas: %d / %d", len(tfs_present), len(tf_list))

    per_species: dict[str, pd.DataFrame] = {}
    for sp in species_levels:
        sp_mask = (sub.obs[cfg.species_key] == sp).values
        if sp_mask.sum() < cfg.min_cells_per_species:
            logger.warning("Skipping %s — only %d cells in lineage", sp, int(sp_mask.sum()))
            continue
        scores = sub.obs.loc[sp_mask, cfg.score_col].astype(float).values
        df = _spearman_block(sub, sp_mask, tfs_present, scores, cfg)
        df["rank"] = df["rho"].abs().rank(method="min", ascending=False)
        per_species[sp] = df

    if not per_species:
        return pd.DataFrame(columns=["tf", "mean_rank", "rank_gap"])

    out = pd.DataFrame({"tf": tfs_present}).set_index("tf")
    for sp, df in per_species.items():
        out[f"{sp}_rho"] = df["rho"].reindex(out.index)
        out[f"{sp}_pval"] = df["pval"].reindex(out.index)
        out[f"{sp}_rank"] = df["rank"].reindex(out.index)

    rank_cols = [c for c in out.columns if c.endswith("_rank")]
    out["mean_rank"] = out[rank_cols].mean(axis=1)
    if len(rank_cols) == 2:
        out["rank_gap"] = out[rank_cols[0]] - out[rank_cols[1]]
    else:
        out["rank_gap"] = 0.0
    return out.dropna(subset=rank_cols, how="all").reset_index()


def _spearman_block(
    sub: ad.AnnData,
    sp_mask: np.ndarray,
    tfs: list[str],
    scores: np.ndarray,
    cfg: RegulatorConfig,
) -> pd.DataFrame:
    rows: list[dict[str, float]] = []
    for tf in tfs:
        x = sub[sp_mask, tf].X
        x = x.toarray().ravel() if hasattr(x, "toarray") else np.asarray(x).ravel()
        if (x > 0).mean() < cfg.min_nonzero_fraction:
            rows.append({"tf": tf, "rho": float("nan"), "pval": float("nan")})
            continue
        try:
            rho, pval = stats.spearmanr(x, scores)
        except (ValueError, FloatingPointError):
            rho, pval = float("nan"), float("nan")
        rows.append({"tf": tf, "rho": float(rho), "pval": float(pval)})
    return pd.DataFrame(rows).set_index("tf")


def split_regulators(
    ranked: pd.DataFrame,
    cap: int = 25,
) -> dict[str, pd.DataFrame]:
    """Return ``{"conserved": …, "human_biased": …, "mouse_biased": …}``.

    Conserved: lowest ``mean_rank`` (high signal in both species), capped.
    Species-biased: largest absolute ``rank_gap`` favouring that species.
    """
    if ranked.empty:
        return {"conserved": ranked, "human_biased": ranked, "mouse_biased": ranked}

    conserved = ranked.sort_values("mean_rank").head(cap).reset_index(drop=True)
    if "human_rank" in ranked.columns and "mouse_rank" in ranked.columns:
        # rank_gap = human_rank - mouse_rank;
        # large positive → high mouse rank but low human → mouse-biased.
        # large negative → high human rank but low mouse → human-biased.
        mouse_biased = (
            ranked.dropna(subset=["rank_gap"])
            .sort_values("rank_gap", ascending=False)
            .head(cap)
            .reset_index(drop=True)
        )
        human_biased = (
            ranked.dropna(subset=["rank_gap"])
            .sort_values("rank_gap", ascending=True)
            .head(cap)
            .reset_index(drop=True)
        )
    else:
        mouse_biased = ranked.head(0)
        human_biased = ranked.head(0)
    return {
        "conserved": conserved,
        "human_biased": human_biased,
        "mouse_biased": mouse_biased,
    }
