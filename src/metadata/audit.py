"""Generate metadata completeness audit report.

Checks all processed h5ad files for harmonized .obs columns and
reports missing values, inconsistencies, and coverage.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import pandas as pd

from src.metadata.annotate import REQUIRED_OBS_COLUMNS

logger = logging.getLogger(__name__)


def audit_metadata(
    processed_dir: Path,
    output_path: Path | None = None,
) -> pd.DataFrame:
    """Audit metadata completeness across all processed h5ad files.

    Parameters
    ----------
    processed_dir : Path
        Directory containing processed h5ad files.
    output_path : Path, optional
        If provided, write the report to this path (CSV or markdown).

    Returns
    -------
    pd.DataFrame
        Audit summary with one row per dataset.
    """
    h5ad_files = sorted(processed_dir.glob("*.h5ad"))
    if not h5ad_files:
        logger.warning("No h5ad files found in %s", processed_dir)
        return pd.DataFrame()

    records = []
    for path in h5ad_files:
        adata = ad.read_h5ad(path, backed="r")
        record = _audit_one(path.stem, adata)
        records.append(record)
        adata.file.close()

    df = pd.DataFrame(records)

    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        if output_path.suffix == ".csv":
            df.to_csv(output_path, index=False)
        else:
            _write_markdown(df, output_path)
        logger.info("Audit report written to %s", output_path)

    return df


def _audit_one(accession: str, adata: ad.AnnData) -> dict:
    """Audit a single dataset's metadata."""
    record = {
        "accession": accession,
        "n_cells": adata.n_obs,
        "n_genes": adata.n_vars,
    }

    obs_cols = set(adata.obs.columns)
    for col in REQUIRED_OBS_COLUMNS:
        present = col in obs_cols
        record[f"has_{col}"] = present
        if present:
            null_frac = adata.obs[col].isna().mean()
            unique_count = adata.obs[col].nunique()
            record[f"{col}_null_pct"] = round(null_frac * 100, 1)
            record[f"{col}_unique"] = unique_count
        else:
            record[f"{col}_null_pct"] = 100.0
            record[f"{col}_unique"] = 0

    return record


def _write_markdown(df: pd.DataFrame, path: Path) -> None:
    """Write audit report as markdown."""
    with open(path, "w") as fh:
        fh.write("# Metadata Audit Report\n\n")
        fh.write(df.to_markdown(index=False))
        fh.write("\n")
