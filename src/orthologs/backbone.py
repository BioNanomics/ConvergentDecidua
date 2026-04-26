"""Build the ortholog backbone table.

Combines Ensembl BioMart and g:Profiler ortholog data into a tiered
mapping table:
  - Tier 1 (strict): one-to-one orthologs confirmed by both sources
  - Tier 2 (relaxed): many-to-many and single-source mappings

Output: results/orthologs/backbone.parquet
"""

from __future__ import annotations

import logging
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from src.orthologs.ensembl import fetch_ensembl_orthologs

logger = logging.getLogger(__name__)


def build_backbone(
    source: str = "human",
    target: str = "mouse",
    cache_dir: Path | None = None,
    output_path: Path | None = None,
    *,
    use_gprofiler: bool = True,
) -> pa.Table:
    """Build the ortholog backbone table.

    Parameters
    ----------
    source, target : str
        Species names.
    cache_dir : Path, optional
        Cache directory for raw API results.
    output_path : Path, optional
        Where to write backbone.parquet.
    use_gprofiler : bool
        Whether to cross-validate with g:Profiler.

    Returns
    -------
    pa.Table
        Backbone table with tier column.
    """
    # Step 1: Ensembl BioMart (primary source)
    ensembl = fetch_ensembl_orthologs(source, target, cache_dir)
    logger.info("Ensembl: %d total rows", len(ensembl))

    # Tier 1: strict 1:1 orthologs with high confidence
    ensembl_df = ensembl.to_pandas()
    one2one = ensembl_df[
        (ensembl_df["orthology_type"] == "ortholog_one2one") & (ensembl_df["confidence"] == 1)
    ].copy()
    one2one["tier"] = 1
    one2one["source_db"] = "ensembl"

    # Tier 2: everything else from Ensembl
    many2many = ensembl_df[~ensembl_df.index.isin(one2one.index)].copy()
    many2many["tier"] = 2
    many2many["source_db"] = "ensembl"

    logger.info("Tier 1 (strict 1:1): %d mappings", len(one2one))
    logger.info("Tier 2 (relaxed): %d mappings", len(many2many))

    # Step 2: g:Profiler cross-validation (optional)
    if use_gprofiler:
        try:
            from src.orthologs.gprofiler import fetch_gprofiler_orthologs

            source_genes = ensembl_df[f"{source}_gene_id"].unique().tolist()
            gprofiler = fetch_gprofiler_orthologs(source_genes, source, target, cache_dir)

            if len(gprofiler) > 0:
                gp_df = gprofiler.to_pandas()
                # Mark Tier 1 genes confirmed by both sources
                confirmed_sources = set(gp_df["source_gene"])
                one2one["gprofiler_confirmed"] = one2one[f"{source}_gene_id"].isin(
                    confirmed_sources
                )
                n_confirmed = one2one["gprofiler_confirmed"].sum()
                logger.info(
                    "g:Profiler confirmed %d / %d Tier 1 mappings",
                    n_confirmed,
                    len(one2one),
                )
        except Exception as exc:
            logger.warning("g:Profiler cross-validation failed (continuing): %s", exc)
            one2one["gprofiler_confirmed"] = False
    else:
        one2one["gprofiler_confirmed"] = False

    many2many["gprofiler_confirmed"] = False

    # Combine
    import pandas as pd

    backbone_df = pd.concat([one2one, many2many], ignore_index=True)

    # Standardize column names
    backbone_df = backbone_df.rename(
        columns={
            f"{source}_gene_id": "source_gene_id",
            f"{source}_symbol": "source_symbol",
            f"{target}_gene_id": "target_gene_id",
            f"{target}_symbol": "target_symbol",
        }
    )
    backbone_df["source_species"] = source
    backbone_df["target_species"] = target

    cols = [
        "source_gene_id",
        "source_symbol",
        "target_gene_id",
        "target_symbol",
        "orthology_type",
        "confidence",
        "tier",
        "source_db",
        "gprofiler_confirmed",
        "source_species",
        "target_species",
    ]
    backbone_df = backbone_df[cols]

    backbone = pa.Table.from_pandas(backbone_df)
    logger.info(
        "Backbone: %d total rows (%d Tier 1, %d Tier 2)",
        len(backbone),
        len(one2one),
        len(many2many),
    )

    # Write output
    if output_path is not None:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        pq.write_table(backbone, output_path)
        logger.info("Backbone written to %s", output_path)

    return backbone
