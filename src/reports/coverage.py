"""Dataset coverage report.

Summarizes which datasets have been successfully processed through
each pipeline stage (fetched, QC'd, integrated, scored).
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from wombat.config import load_config

logger = logging.getLogger(__name__)


def generate_coverage_report(results_dir: Path, output_path: Path) -> pd.DataFrame:
    """Generate dataset coverage report.

    Parameters
    ----------
    results_dir : Path
        Path to the results directory.
    output_path : Path
        Where to write the report.

    Returns
    -------
    pd.DataFrame
        Coverage matrix.
    """
    datasets = load_config("datasets")
    records = []

    for ds in datasets:
        acc = ds["accession"]
        record = {
            "accession": acc,
            "species": ds["species"],
            "assay": ds["assay"],
            "fetched": (results_dir / "processed" / f"{acc}.h5ad").exists(),
            "qc_passed": (results_dir / "qc" / f"{acc}.h5ad").exists(),
            "integrated": False,  # checked below
            "scored": False,  # checked below
        }
        records.append(record)

    # Check if integrated/scored data exists
    integrated_path = results_dir / "integrated" / "stromal_harmony.h5ad"
    scored_path = results_dir / "scored" / "stromal_scored.h5ad"

    df = pd.DataFrame(records)

    if integrated_path.exists():
        df.loc[df["assay"].str.contains("scRNA|snRNA", case=False, na=False), "integrated"] = True

    if scored_path.exists():
        df.loc[df["assay"].str.contains("scRNA|snRNA", case=False, na=False), "scored"] = True

    # Write report
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        fh.write("# Dataset Coverage Report\n\n")
        fh.write(df.to_markdown(index=False))
        fh.write("\n\n## Summary\n\n")
        for col in ["fetched", "qc_passed", "integrated", "scored"]:
            n = df[col].sum()
            fh.write(f"- **{col}**: {n}/{len(df)} datasets\n")

    logger.info("Coverage report → %s", output_path)
    return df
