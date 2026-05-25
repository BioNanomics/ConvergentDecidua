"""Cross-species stromal integration.

Maps gene space to the ortholog backbone, then runs Harmony (default)
or scVI to produce a joint embedding of human and mouse stromal cells.
"""

from __future__ import annotations

import logging
from pathlib import Path

import anndata as ad
import pandas as pd
import pyarrow.parquet as pq
import scanpy as sc

logger = logging.getLogger(__name__)


def integrate_stromal(
    adata_list: list[ad.AnnData],
    backbone_path: Path,
    method: str = "harmony",
    n_hvg: int = 3000,
) -> ad.AnnData:
    """Integrate stromal datasets across species.

    Parameters
    ----------
    adata_list : list[ad.AnnData]
        List of per-dataset stromal AnnData objects (already QC'd and
        annotated with ``species``, ``dataset``, ``cell_type`` in .obs).
    backbone_path : Path
        Path to ``results/orthologs/backbone.parquet``.
    method : str
        Integration method: ``"harmony"`` (default) or ``"scvi"``.
    n_hvg : int
        Number of highly variable genes to select.

    Returns
    -------
    ad.AnnData
        Integrated AnnData with joint embedding in ``.obsm["X_integrated"]``.
    """
    # Step 1: Map all datasets to common gene space via backbone
    backbone = pq.read_table(backbone_path).to_pandas()
    tier1 = backbone[backbone["tier"] == 1]

    mapped = []
    for adata in adata_list:
        species = adata.obs["species"].iloc[0]
        if species == "human":
            mapped.append(_subset_to_backbone_genes(adata, tier1, "source"))
        else:
            mapped.append(_remap_mouse_genes(adata, tier1))

    # Step 2: Concatenate
    combined = ad.concat(
        mapped, join="inner", label="batch", keys=[f"{a.obs['dataset'].iloc[0]}" for a in mapped]
    )
    logger.info("Combined: %d cells × %d genes", combined.n_obs, combined.n_vars)

    # Step 3: Re-normalize and HVG
    if "counts" in combined.layers:
        combined.X = combined.layers["counts"].copy()
    sc.pp.normalize_total(combined, target_sum=1e4)
    sc.pp.log1p(combined)
    sc.pp.highly_variable_genes(combined, n_top_genes=min(n_hvg, combined.n_vars))
    combined = combined[:, combined.var["highly_variable"]].copy()
    sc.pp.scale(combined, max_value=10)
    sc.tl.pca(combined, n_comps=min(50, combined.n_obs - 1, combined.n_vars - 1))

    # Step 4: Integrate
    if method == "harmony":
        _run_harmony(combined)
    elif method == "scvi":
        _run_scvi(combined)
    else:
        msg = f"Unknown integration method: {method}"
        raise ValueError(msg)

    # Step 5: Neighbors and UMAP on corrected embedding
    rep_key = "X_pca_harmony" if method == "harmony" else "X_scvi"
    sc.pp.neighbors(combined, use_rep=rep_key)
    sc.tl.umap(combined)
    combined.obsm["X_integrated"] = combined.obsm[rep_key]

    logger.info("Integration complete (%s): %d cells", method, combined.n_obs)
    return combined


def _subset_to_backbone_genes(
    adata: ad.AnnData,
    tier1: pd.DataFrame,
    col_prefix: str,
) -> ad.AnnData:
    """Subset human AnnData to genes in the backbone."""
    backbone_genes = set(tier1[f"{col_prefix}_symbol"])
    present = [g for g in adata.var_names if g in backbone_genes]
    return adata[:, present].copy()


def _remap_mouse_genes(adata: ad.AnnData, tier1: pd.DataFrame) -> ad.AnnData:
    """Remap mouse gene symbols to human orthologs via backbone."""
    mouse_to_human = dict(zip(tier1["target_symbol"], tier1["source_symbol"], strict=False))

    # Filter to genes with orthologs
    has_ortholog = [g for g in adata.var_names if g in mouse_to_human]
    adata = adata[:, has_ortholog].copy()

    # Rename to human symbols
    adata.var_names = pd.Index([mouse_to_human[g] for g in adata.var_names])

    # Handle duplicates (take the first)
    if adata.var_names.duplicated().any():
        adata = adata[:, ~adata.var_names.duplicated()].copy()

    return adata


def _run_harmony(adata: ad.AnnData) -> None:
    """Run Harmony integration on PCA coordinates."""
    try:
        import harmonypy

        # Use 'dataset' as batch key if multiple datasets, 'species' if multi-species
        n_species = adata.obs["species"].nunique()
        batch_key = "species" if n_species > 1 else "dataset"
        logger.info("Harmony batch key: %s (%d unique)", batch_key, adata.obs[batch_key].nunique())

        ho = harmonypy.run_harmony(
            adata.obsm["X_pca"],
            adata.obs,
            batch_key,
            max_iter_harmony=20,
        )
        corrected = ho.Z_corr.T
        # Ensure shape matches (n_cells, n_pcs)
        if corrected.shape[0] != adata.n_obs:
            corrected = ho.Z_corr if ho.Z_corr.shape[0] == adata.n_obs else corrected.T
        adata.obsm["X_pca_harmony"] = corrected
        logger.info("Harmony integration complete")
    except ImportError:
        logger.warning("harmonypy not installed — falling back to uncorrected PCA")
        adata.obsm["X_pca_harmony"] = adata.obsm["X_pca"]


def _run_scvi(adata: ad.AnnData) -> None:
    """Run scVI integration (requires scvi-tools)."""
    try:
        import scvi as scvi_tools

        scvi_tools.model.SCVI.setup_anndata(adata, batch_key="species", layer="counts")
        model = scvi_tools.model.SCVI(adata, n_latent=30)
        model.train(max_epochs=100, early_stopping=True)
        adata.obsm["X_scvi"] = model.get_latent_representation()
        logger.info("scVI integration complete")
    except ImportError as exc:
        msg = "scvi-tools not installed — install with: pip install scvi-tools"
        raise ImportError(msg) from exc
