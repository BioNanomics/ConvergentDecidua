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

    df = pd.DataFrame(records)

    # Honestly detect which datasets actually contributed cells to the
    # integrated / scored objects by reading .obs.dataset, rather than
    # blanket-marking every scRNA dataset whenever the file exists.
    integrated_path = results_dir / "integrated" / "stromal_harmony.h5ad"
    scored_path = results_dir / "scored" / "stromal_scored.h5ad"

    integrated_accessions = _accessions_in_h5ad(integrated_path)
    if integrated_accessions is not None:
        df["integrated"] = df["accession"].isin(integrated_accessions)

    scored_accessions = _accessions_in_h5ad(scored_path)
    if scored_accessions is not None:
        df["scored"] = df["accession"].isin(scored_accessions)

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


def _accessions_in_h5ad(path: Path) -> set[str] | None:
    """Return the set of dataset accessions present in an h5ad's ``.obs.dataset``.

    Returns ``None`` when the file does not exist or cannot be read; returns an
    empty set when the file exists but lacks the ``dataset`` column.
    """
    if not path.exists():
        return None
    try:
        import anndata as ad

        a = ad.read_h5ad(path, backed="r")
        try:
            if "dataset" not in a.obs.columns:
                return set()
            return set(a.obs["dataset"].astype(str).unique())
        finally:
            if getattr(a, "file", None) is not None:
                a.file.close()
    except Exception as exc:  # pragma: no cover - defensive
        logger.warning("Could not inspect %s for coverage detection: %s", path, exc)
        return None
