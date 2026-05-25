"""Marker-based cell-type annotation.

Scores each cell against the cell-type marker gene sets from
configs/markers.yaml and assigns the best-matching label.
"""

from __future__ import annotations

import logging

import anndata as ad
import scanpy as sc

from wombat.config import load_config

logger = logging.getLogger(__name__)


def annotate_cell_types(
    adata: ad.AnnData,
    species: str = "human",
    backbone_path: str | None = None,
) -> ad.AnnData:
    """Score cells and assign cell-type labels.

    Parameters
    ----------
    adata : ad.AnnData
        QC'd, log-normalized AnnData.
    species : str
        Species name. If not human, gene sets are mapped via backbone.
    backbone_path : str, optional
        Path to ortholog backbone parquet (required for non-human species).

    Returns
    -------
    ad.AnnData
        AnnData with ``cell_type`` and ``cell_type_score`` in ``.obs``.
    """
    markers_cfg = load_config("markers")
    cell_type_markers = markers_cfg["cell_type_markers"]

    gene_sets = _prepare_gene_sets(cell_type_markers, species, backbone_path, adata.var_names)

    if not gene_sets:
        logger.warning("No usable gene sets after filtering — skipping annotation")
        adata.obs["cell_type"] = "unknown"
        adata.obs["cell_type_score"] = 0.0
        return adata

    # Score each cell type
    score_cols = []
    for ct, genes in gene_sets.items():
        col = f"score_{ct}"
        sc.tl.score_genes(adata, gene_list=genes, score_name=col)
        score_cols.append((ct, col))

    # Assign best label
    import pandas as pd

    score_df = pd.DataFrame(
        {ct: adata.obs[col] for ct, col in score_cols},
        index=adata.obs.index,
    )
    adata.obs["cell_type"] = score_df.idxmax(axis=1)
    adata.obs["cell_type_score"] = score_df.max(axis=1)

    # Log distribution
    counts = adata.obs["cell_type"].value_counts()
    logger.info("Cell-type annotation: %d types assigned", len(counts))
    for ct, n in counts.items():
        logger.debug("  %s: %d cells", ct, n)

    return adata


def _prepare_gene_sets(
    markers: dict[str, list[str]],
    species: str,
    backbone_path: str | None,
    var_names: pd.Index,
) -> dict[str, list[str]]:
    """Prepare gene sets, mapping through backbone if needed."""
    if species != "human" and backbone_path:
        symbol_map = _load_symbol_map(backbone_path, direction="human_to_mouse")
    else:
        symbol_map = None

    var_set = set(var_names)
    result = {}

    from src.scoring.gene_sets import apply_species_overrides

    for ct, genes in markers.items():
        if not genes:
            continue

        mapped = [symbol_map.get(g, g) for g in genes] if symbol_map else list(genes)

        # Augment with per-species overrides for genes the backbone cannot map
        # (e.g. mouse decidual-prolactin family; PRL is Tier 2 only).
        mapped = apply_species_overrides(ct, mapped, species, "cell_type_markers")

        # Filter to genes present in the dataset
        present = [g for g in mapped if g in var_set]
        if len(present) >= 2:
            result[ct] = present
        else:
            logger.debug("Skipping %s: only %d/%d genes found", ct, len(present), len(mapped))

    return result


def _load_symbol_map(backbone_path: str, direction: str = "human_to_mouse") -> dict[str, str]:
    """Load symbol mapping from backbone parquet."""
    import pyarrow.parquet as pq

    table = pq.read_table(backbone_path)
    df = table.to_pandas()

    # Tier 1 only for annotation
    tier1 = df[df["tier"] == 1]

    if direction == "human_to_mouse":
        return dict(zip(tier1["source_symbol"], tier1["target_symbol"], strict=False))
    return dict(zip(tier1["target_symbol"], tier1["source_symbol"], strict=False))


import pandas as pd  # noqa: E402  — needed for pd.Index type in _prepare_gene_sets
