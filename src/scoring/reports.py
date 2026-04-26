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
