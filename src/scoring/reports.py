"""Generate scoring report figures.

Produces heatmaps and violin plots of decidualization scores
by cell state and species.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import pandas as pd

logger = logging.getLogger(__name__)


def generate_score_report(
    adata: ad.AnnData,
    score_columns: list[str],
    output_dir: Path,
) -> Path:
    """Generate scoring report with heatmap and violin plots.

    Parameters
    ----------
    adata : ad.AnnData
        Scored, integrated AnnData.
    score_columns : list[str]
        Column names in .obs that contain scores.
    output_dir : Path
        Directory to write report files.

    Returns
    -------
    Path
        Path to the generated report directory.
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Summary table: mean scores by cell_type × species
    summary = _build_summary(adata, score_columns)
    summary.to_csv(output_dir / "score_summary.csv")
    logger.info("Score summary → %s", output_dir / "score_summary.csv")

    # Heatmap
    _plot_heatmap(summary, output_dir / "score_heatmap.png")

    # Violin plots
    _plot_violins(adata, score_columns, output_dir)

    # Markdown report
    _write_markdown_report(summary, score_columns, output_dir / "scoring_report.md")

    logger.info("Scoring report written to %s", output_dir)
    return output_dir


def _build_summary(adata: ad.AnnData, score_columns: list[str]) -> pd.DataFrame:
    """Build mean score summary grouped by cell_type and species."""
    groupby_cols = []
    if "cell_type" in adata.obs.columns:
        groupby_cols.append("cell_type")
    if "species" in adata.obs.columns:
        groupby_cols.append("species")

    if not groupby_cols:
        groupby_cols = ["dataset"]

    present_scores = [c for c in score_columns if c in adata.obs.columns]
    return adata.obs.groupby(groupby_cols)[present_scores].mean()


def _plot_heatmap(summary: pd.DataFrame, output_path: Path) -> None:
    """Plot score heatmap."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns

        fig, ax = plt.subplots(figsize=(12, max(6, len(summary) * 0.4)))
        sns.heatmap(summary, annot=True, fmt=".2f", cmap="RdYlBu_r", ax=ax)
        ax.set_title("Decidualization Scores by Cell Type × Species")
        plt.tight_layout()
        fig.savefig(output_path, dpi=150)
        plt.close(fig)
        logger.info("Heatmap → %s", output_path)
    except ImportError:
        logger.warning("matplotlib/seaborn not available — skipping heatmap")


def _plot_violins(adata: ad.AnnData, score_columns: list[str], output_dir: Path) -> None:
    """Plot violin plots for each score module."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns

        groupby = "cell_type" if "cell_type" in adata.obs.columns else "dataset"

        for col in score_columns:
            if col not in adata.obs.columns:
                continue

            fig, ax = plt.subplots(figsize=(10, 5))
            plot_df = adata.obs[[groupby, col]].copy()
            if "species" in adata.obs.columns:
                plot_df["species"] = adata.obs["species"]
                sns.violinplot(data=plot_df, x=groupby, y=col, hue="species", ax=ax, split=True)
            else:
                sns.violinplot(data=plot_df, x=groupby, y=col, ax=ax)
            ax.set_title(col)
            plt.xticks(rotation=45, ha="right")
            plt.tight_layout()
            fig.savefig(output_dir / f"violin_{col}.png", dpi=150)
            plt.close(fig)

        logger.info("Violin plots → %s", output_dir)
    except ImportError:
        logger.warning("matplotlib/seaborn not available — skipping violin plots")


def _write_markdown_report(
    summary: pd.DataFrame,
    score_columns: list[str],
    output_path: Path,
) -> None:
    """Write scoring report as markdown."""
    with open(output_path, "w") as fh:
        fh.write("# Decidualization Scoring Report\n\n")
        fh.write("## Score Summary\n\n")
        fh.write(summary.to_markdown())
        fh.write("\n\n## Modules Scored\n\n")
        for col in score_columns:
            fh.write(f"- {col}\n")
        fh.write("\n## Figures\n\n")
        fh.write("- `score_heatmap.png` — Mean scores by cell type × species\n")
        for col in score_columns:
            fh.write(f"- `violin_{col}.png` — Distribution by cell type\n")
    logger.info("Markdown report → %s", output_path)


def generate_bulk_score_report(
    scored: dict[str, ad.AnnData],
    monotonicity_tables: dict[str, pd.DataFrame],
    output_dir: Path,
) -> Path:
    """Write a combined bulk-scoring report: per-dataset monotonicity tables
    and a per-module ``decidual_score`` vs. time line plot for each dataset.

    Parameters
    ----------
    scored
        Mapping ``accession → scored AnnData`` (with numeric ``time`` in obs
        and one column per scoring module).
    monotonicity_tables
        Mapping ``accession → DataFrame`` returned by
        ``src.scoring.bulk.monotonicity``.
    output_dir
        Directory to write the report + plots into.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    md_path = output_dir / "bulk_scoring_report.md"

    _plot_bulk_score_vs_time(scored, output_dir)

    with open(md_path, "w") as fh:
        fh.write("# Bulk RNA-seq scoring — monotonicity\n\n")
        fh.write(
            "For each bulk dataset, Spearman rank correlation between the\n"
            "experimental ``time`` axis (parsed from sample labels) and\n"
            "each scoring module. ``monotonic`` = |rho| ≥ 0.7 and pval <\n"
            "0.05 (Q3.1 floor; small-n bulk has limited statistical power\n"
            "so rho carries the signal). A monotone ``decidual_score`` is\n"
            "the Q3.1 acceptance criterion.\n\n"
        )
        for acc, table in monotonicity_tables.items():
            adata = scored[acc]
            fh.write(f"## {acc}\n\n")
            fh.write(
                f"- samples: **{adata.n_obs}**, "
                f"genes: **{adata.n_vars}**, "
                f"species: **{adata.obs['species'].iloc[0]}**\n"
            )
            fh.write(
                f"- time axis: "
                f"`{dict(zip(adata.obs_names, adata.obs['time'].astype(float), strict=False))}`\n\n"
            )
            fh.write(table.reset_index().to_markdown(index=False))
            fh.write(f"\n\n![{acc} decidual_score vs time](decidual_score_vs_time_{acc}.png)\n\n")
    logger.info("Bulk scoring report → %s", md_path)
    return md_path


def _plot_bulk_score_vs_time(scored: dict[str, ad.AnnData], output_dir: Path) -> None:
    """Plot every module vs. time for each bulk dataset.

    Highlights ``decidual_score`` (the Q3.1 acceptance signal) in colour;
    other modules drawn in light grey for context.
    """
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        logger.warning("matplotlib not available — skipping bulk score plots")
        return

    for acc, adata in scored.items():
        if "time" not in adata.obs.columns:
            continue
        time = adata.obs["time"].astype(float)
        order = time.argsort().values
        time = time.iloc[order]
        score_cols = [c for c in adata.obs.columns if c.endswith("_score")]
        if not score_cols:
            continue

        fig, ax = plt.subplots(figsize=(8, 5))
        for col in score_cols:
            y = adata.obs[col].astype(float).iloc[order]
            if col == "decidual_score":
                ax.plot(time, y, marker="o", linewidth=2, label=col, color="C3", zorder=5)
            else:
                ax.plot(time, y, marker=".", linewidth=1, alpha=0.5, color="grey", label=col)
        ax.set_xlabel("time (parsed from sample labels)")
        ax.set_ylabel("module score (z-scored mean)")
        ax.set_title(f"{acc} — module scores vs time")
        ax.axhline(0, color="black", linewidth=0.5, linestyle=":")
        ax.legend(loc="best", fontsize="x-small", ncol=2)
        fig.tight_layout()
        out = output_dir / f"decidual_score_vs_time_{acc}.png"
        fig.savefig(out, dpi=150)
        plt.close(fig)
        logger.info("Plot → %s", out)
