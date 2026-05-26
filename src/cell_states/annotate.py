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
    cell_type_lineages = markers_cfg.get("cell_type_lineages", {})

    gene_sets = _prepare_gene_sets(cell_type_markers, species, backbone_path, adata.var_names)

    if not gene_sets:
        logger.warning("No usable gene sets after filtering — skipping annotation")
        adata.obs["cell_type"] = "unknown"
        adata.obs["cell_type_score"] = 0.0
        adata.obs["lineage"] = "unknown"
        return adata

    # Score each cell type
    score_cols = []
    for ct, genes in gene_sets.items():
        col = f"score_{ct}"
        sc.tl.score_genes(adata, gene_list=genes, score_name=col)
        score_cols.append((ct, col))

    import pandas as pd

    score_df = pd.DataFrame(
        {ct: adata.obs[col] for ct, col in score_cols},
        index=adata.obs.index,
    )

    if cell_type_lineages:
        # Hierarchical assignment: lineage first (max over constituent
        # cell-type scores), then fine-grained cell_type within winning
        # lineage. Avoids vote-splitting across stromal sub-types.
        adata.obs["cell_type"], adata.obs["cell_type_score"], adata.obs["lineage"] = (
            _hierarchical_assign(score_df, cell_type_lineages)
        )
    else:
        # Legacy flat idxmax (no lineages declared in markers.yaml)
        adata.obs["cell_type"] = score_df.idxmax(axis=1)
        adata.obs["cell_type_score"] = score_df.max(axis=1)
        adata.obs["lineage"] = "unknown"

    # Log distribution
    counts = adata.obs["cell_type"].value_counts()
    logger.info("Cell-type annotation: %d types assigned", len(counts))
    for ct, n in counts.items():
        logger.debug("  %s: %d cells", ct, n)

    return adata


def _hierarchical_assign(
    score_df,
    lineages: dict[str, list[str]],
):
    """Two-pass assignment: lineage by max sub-score, then cell_type within.

    Returns a tuple of (cell_type, cell_type_score, lineage) Series, all
    indexed like ``score_df``.
    """
    import pandas as pd

    # Build lineage→score matrix: each lineage's score per cell is the
    # max of its constituent cell-type scores that are actually present
    # in ``score_df``.
    lineage_scores = {}
    members: dict[str, list[str]] = {}
    for lineage_name, ct_list in lineages.items():
        present = [ct for ct in ct_list if ct in score_df.columns]
        if not present:
            continue
        members[lineage_name] = present
        lineage_scores[lineage_name] = score_df[present].max(axis=1)

    if not lineage_scores:
        # No declared lineage members were scored — fall back to flat idxmax.
        return (
            score_df.idxmax(axis=1),
            score_df.max(axis=1),
            pd.Series("unknown", index=score_df.index),
        )

    lineage_df = pd.DataFrame(lineage_scores, index=score_df.index)
    lineage = lineage_df.idxmax(axis=1)

    # Within winning lineage, pick the fine-grained cell_type.
    cell_type = pd.Series(index=score_df.index, dtype=object)
    for lineage_name, present in members.items():
        mask = lineage == lineage_name
        if not mask.any():
            continue
        sub = score_df.loc[mask, present]
        cell_type.loc[mask] = sub.idxmax(axis=1)

    cell_type_score = pd.Series(
        [score_df.loc[idx, ct] for idx, ct in cell_type.items()],
        index=score_df.index,
    )

    return cell_type, cell_type_score, lineage


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
