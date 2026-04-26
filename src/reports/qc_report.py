"""QC summary report.

Aggregates QC metrics across all processed datasets into a single report.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def generate_qc_report(results_dir: Path, output_path: Path) -> pd.DataFrame:
    """Generate QC summary report.

    Parameters
    ----------
    results_dir : Path
        Path to the results directory.
    output_path : Path
        Where to write the report.

    Returns
    -------
    pd.DataFrame
        QC summary with one row per dataset.
    """
    qc_dir = results_dir / "qc"
    h5ad_files = sorted(qc_dir.glob("*.h5ad")) if qc_dir.exists() else []

    if not h5ad_files:
        logger.warning("No QC'd h5ad files found")
        return pd.DataFrame()

    import anndata as ad

    records = []
    for path in h5ad_files:
        adata = ad.read_h5ad(path, backed="r")
        record = {
            "accession": path.stem,
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
        }

        if "species" in adata.obs.columns:
            record["species"] = adata.obs["species"].iloc[0]
        if "n_genes_by_counts" in adata.obs.columns:
            record["median_genes"] = int(adata.obs["n_genes_by_counts"].median())
        if "pct_counts_mt" in adata.obs.columns:
            record["median_pct_mito"] = round(adata.obs["pct_counts_mt"].median(), 1)
        if "qc_passed" in adata.obs.columns:
            record["all_qc_passed"] = bool(adata.obs["qc_passed"].all())

        records.append(record)
        adata.file.close()

    df = pd.DataFrame(records)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        fh.write("# QC Summary Report\n\n")
        fh.write(df.to_markdown(index=False))
        fh.write("\n")

    logger.info("QC report → %s", output_path)
    return df
