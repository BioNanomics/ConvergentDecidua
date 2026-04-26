"""Extract stromal-lineage cell subsets.

Subsets QC'd AnnData to the 4 stromal subtypes defined in the
cell-state ontology: stromal_fibroblast, pre_decidual_stromal,
decidual_stromal, senescent_decidual.
"""

from __future__ import annotations

import logging

import anndata as ad

logger = logging.getLogger(__name__)

STROMAL_TYPES = [
    "stromal_fibroblast",
    "pre_decidual_stromal",
    "decidual_stromal",
    "senescent_decidual",
]


def subset_stromal(adata: ad.AnnData) -> ad.AnnData:
    """Extract stromal cells from an annotated AnnData.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData with ``cell_type`` column in ``.obs``.

    Returns
    -------
    ad.AnnData
        Subset containing only stromal cells.
    """
    if "cell_type" not in adata.obs.columns:
        msg = "AnnData missing 'cell_type' column — run annotate first"
        raise ValueError(msg)

    mask = adata.obs["cell_type"].isin(STROMAL_TYPES)
    stromal = adata[mask].copy()
    logger.info(
        "Stromal subset: %d / %d cells (%d types)",
        stromal.n_obs,
        adata.n_obs,
        stromal.obs["cell_type"].nunique(),
    )
    return stromal
