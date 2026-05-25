"""Generic gene-set scoring engine.

Wraps scanpy.tl.score_genes with species-aware gene-set mapping
through the ortholog backbone.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import pyarrow.parquet as pq
import scanpy as sc

logger = logging.getLogger(__name__)


def score_gene_set(
    adata: ad.AnnData,
    genes: list[str],
    score_name: str,
    species: str = "human",
    backbone_path: Path | None = None,
) -> ad.AnnData:
    """Score cells for a single gene set.

    Parameters
    ----------
    adata : ad.AnnData
        Log-normalized AnnData.
    genes : list[str]
        Human gene symbols.
    score_name : str
        Column name for the score in ``.obs``.
    species : str
        Species of the dataset.
    backbone_path : Path, optional
        Path to ortholog backbone (required for non-human).

    Returns
    -------
    ad.AnnData
        AnnData with score added to ``.obs[score_name]``.
    """
    mapped = _map_genes(genes, species, backbone_path, adata.var_names, set_name=score_name)

    if len(mapped) < 2:
        logger.warning(
            "Score '%s': only %d genes found — setting to 0",
            score_name,
            len(mapped),
        )
        adata.obs[score_name] = 0.0
        return adata

    sc.tl.score_genes(adata, gene_list=mapped, score_name=score_name)
    logger.info(
        "Scored '%s': %d/%d genes used, mean=%.3f",
        score_name,
        len(mapped),
        len(genes),
        adata.obs[score_name].mean(),
    )
    return adata


def score_all_modules(
    adata: ad.AnnData,
    gene_sets: dict[str, list[str]],
    species: str = "human",
    backbone_path: Path | None = None,
) -> ad.AnnData:
    """Score cells for all gene-set modules.

    Parameters
    ----------
    adata : ad.AnnData
        Log-normalized AnnData.
    gene_sets : dict
        Mapping of score name → gene list (human symbols).
    species : str
        Species of the dataset.
    backbone_path : Path, optional
        Path to ortholog backbone.

    Returns
    -------
    ad.AnnData
        AnnData with all score columns added.
    """
    for name, genes in gene_sets.items():
        adata = score_gene_set(adata, genes, name, species, backbone_path)
    return adata


def _map_genes(
    genes: list[str],
    species: str,
    backbone_path: Path | None,
    var_names,
    set_name: str | None = None,
) -> list[str]:
    """Map gene symbols to the target species and filter to present genes.

    When ``set_name`` is provided and ``species != 'human'``, per-species
    overrides from ``markers.yaml::species_overrides`` are merged in before
    filtering against ``var_names``.
    """
    var_set = set(var_names)

    if species != "human" and backbone_path:
        table = pq.read_table(backbone_path)
        df = table.to_pandas()
        tier1 = df[df["tier"] == 1]
        h2m = dict(zip(tier1["source_symbol"], tier1["target_symbol"]))
        mapped = [h2m.get(g, g) for g in genes]
    else:
        mapped = list(genes)

    if set_name and species != "human":
        from src.scoring.gene_sets import apply_species_overrides

        mapped = apply_species_overrides(set_name, mapped, species, "score_gene_sets")

    return [g for g in mapped if g in var_set]
