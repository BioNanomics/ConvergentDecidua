"""Apply harmonized metadata to AnnData objects.

Reads a processed h5ad, adds standardized .obs columns, and writes it back.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad

from src.metadata.harmonize import (
    normalize_assay,
    normalize_cycle_stage,
    normalize_donor,
    normalize_sample,
    normalize_species,
)

logger = logging.getLogger(__name__)

# Required harmonized columns in .obs
REQUIRED_OBS_COLUMNS = ["species", "assay", "dataset", "cycle_stage", "donor", "sample"]


def annotate_h5ad(
    h5ad_path: Path,
    dataset_meta: dict,
    *,
    overwrite: bool = True,
) -> ad.AnnData:
    """Add harmonized metadata columns to an h5ad file.

    Parameters
    ----------
    h5ad_path : Path
        Path to the h5ad file to annotate.
    dataset_meta : dict
        Entry from datasets.yaml for this dataset.
    overwrite : bool
        If True, write the annotated file back to the same path.

    Returns
    -------
    ad.AnnData
        Annotated AnnData object.
    """
    logger.info("Annotating %s", h5ad_path)
    adata = ad.read_h5ad(h5ad_path)
    acc = dataset_meta["accession"]

    # Species — always overwrite from config (authoritative)
    adata.obs["species"] = normalize_species(dataset_meta["species"])

    # Assay
    adata.obs["assay"] = normalize_assay(dataset_meta["assay"])

    # Dataset accession
    adata.obs["dataset"] = acc

    # Cycle stage — try to find in existing .obs, else use condition
    if "cycle_stage" in adata.obs.columns:
        adata.obs["cycle_stage"] = adata.obs["cycle_stage"].map(normalize_cycle_stage)
    elif "phase" in adata.obs.columns:
        adata.obs["cycle_stage"] = adata.obs["phase"].map(normalize_cycle_stage)
    else:
        adata.obs["cycle_stage"] = normalize_cycle_stage(dataset_meta.get("condition", "unknown"))

    # Donor — try to find in existing .obs
    if "donor" in adata.obs.columns:
        adata.obs["donor"] = adata.obs["donor"].apply(lambda x: normalize_donor(str(x), acc))
    elif "patient" in adata.obs.columns:
        adata.obs["donor"] = adata.obs["patient"].apply(lambda x: normalize_donor(str(x), acc))
    else:
        adata.obs["donor"] = normalize_donor("unknown", acc)

    # Sample — try to find in existing .obs
    if "sample" in adata.obs.columns:
        adata.obs["sample"] = adata.obs["sample"].apply(lambda x: normalize_sample(str(x), acc))
    elif "sample_id" in adata.obs.columns:
        adata.obs["sample"] = adata.obs["sample_id"].apply(lambda x: normalize_sample(str(x), acc))
    else:
        adata.obs["sample"] = normalize_sample("unknown", acc)

    if overwrite:
        adata.write_h5ad(h5ad_path)
        logger.info("Updated %s with harmonized metadata", h5ad_path)

    return adata
