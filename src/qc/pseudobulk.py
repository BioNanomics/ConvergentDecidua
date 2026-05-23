"""Generate pseudobulk expression profiles from single-cell data.

Aggregates scRNA-seq cells by sample/donor/cell_type to produce
pseudobulk AnnData objects for downstream bulk-style analyses.
"""

from __future__ import annotations

import logging
from collections import defaultdict

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

logger = logging.getLogger(__name__)

DEFAULT_MIN_CELLS = 20


def pseudobulk(
    adata: ad.AnnData,
    groupby: list[str] | str = "sample",
    layer: str | None = "counts",
    min_cells: int = DEFAULT_MIN_CELLS,
) -> ad.AnnData:
    """Aggregate single-cell data to pseudobulk.

    Parameters
    ----------
    adata : ad.AnnData
        Single-cell AnnData (should have raw counts in layer).
    groupby : str or list[str]
        Column(s) in .obs to group by. E.g. ``"sample"`` or
        ``["sample", "cell_type"]``.
    layer : str, optional
        Layer to aggregate. ``None`` uses .X.
    min_cells : int
        Minimum cells per group to include.

    Returns
    -------
    ad.AnnData
        Pseudobulk AnnData (groups × genes) with summed counts.
    """
    if isinstance(groupby, str):
        groupby = [groupby]

    X = adata.layers[layer] if layer and layer in adata.layers else adata.X
    if sp.issparse(X):
        X = X.toarray()

    # Create group labels
    group_labels = adata.obs[groupby].apply(lambda row: "|".join(row.astype(str)), axis=1)

    # Aggregate
    groups: dict[str, list[int]] = defaultdict(list)
    for i, label in enumerate(group_labels):
        groups[label].append(i)

    kept_labels = []
    summed_rows = []
    obs_records = []

    for label, indices in sorted(groups.items()):
        if len(indices) < min_cells:
            logger.debug(
                "Skipping group %s: only %d cells (min=%d)", label, len(indices), min_cells
            )
            continue
        kept_labels.append(label)
        summed_rows.append(X[indices].sum(axis=0))
        # Record group metadata
        parts = label.split("|")
        record = dict(zip(groupby, parts))
        record["n_cells"] = len(indices)
        obs_records.append(record)

    if not summed_rows:
        msg = f"No groups with >= {min_cells} cells"
        raise ValueError(msg)

    mat = np.vstack(summed_rows)
    obs_df = pd.DataFrame(obs_records, index=kept_labels)
    var_df = adata.var.copy()

    pb = ad.AnnData(X=mat, obs=obs_df, var=var_df)
    logger.info(
        "Pseudobulk: %d groups × %d genes (from %d cells)", pb.n_obs, pb.n_vars, adata.n_obs
    )
    return pb
