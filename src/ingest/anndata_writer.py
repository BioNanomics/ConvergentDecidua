"""Convert downloaded data files to standardized AnnData h5ad format.

Supports multiple input formats: h5ad, 10X MTX, CSV/TSV count matrices.
Attaches metadata from datasets.yaml to .obs and .uns.
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.io
import scipy.sparse

logger = logging.getLogger(__name__)


def to_anndata(
    raw_dir: Path,
    dataset_meta: dict,
    output_path: Path,
) -> ad.AnnData:
    """Convert downloaded files in raw_dir to a standardized h5ad.

    Parameters
    ----------
    raw_dir : Path
        Directory containing downloaded files for one dataset.
    dataset_meta : dict
        Entry from datasets.yaml for this dataset.
    output_path : Path
        Where to write the output h5ad.

    Returns
    -------
    ad.AnnData
        The loaded and annotated AnnData object.
    """
    adata = _load_from_dir(raw_dir)

    # Attach dataset metadata to .uns
    adata.uns["dataset"] = {
        "accession": dataset_meta["accession"],
        "species": dataset_meta["species"],
        "assay": dataset_meta["assay"],
        "tissue": dataset_meta.get("tissue", ""),
        "condition": dataset_meta.get("condition", ""),
    }

    # Ensure .obs has species column
    if "species" not in adata.obs.columns:
        adata.obs["species"] = dataset_meta["species"]
    if "dataset" not in adata.obs.columns:
        adata.obs["dataset"] = dataset_meta["accession"]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_path)
    logger.info("Wrote %d cells × %d genes → %s", adata.n_obs, adata.n_vars, output_path)
    return adata


def _load_from_dir(raw_dir: Path) -> ad.AnnData:
    """Auto-detect format and load AnnData from a directory."""
    # Priority 1: existing h5ad
    h5ad_files = list(raw_dir.glob("*.h5ad")) + list(raw_dir.glob("**/*.h5ad"))
    if h5ad_files:
        logger.info("Loading h5ad: %s", h5ad_files[0])
        return ad.read_h5ad(h5ad_files[0])

    # Priority 2: 10X MTX format (matrix.mtx + barcodes + features/genes)
    mtx_files = list(raw_dir.glob("**/matrix.mtx*"))
    if mtx_files:
        return _load_10x_mtx(mtx_files[0].parent)

    # Priority 3: h5 (10X HDF5)
    h5_files = list(raw_dir.glob("*.h5")) + list(raw_dir.glob("**/*.h5"))
    if h5_files:
        logger.info("Loading 10X h5: %s", h5_files[0])
        return ad.read_10x_h5(h5_files[0])

    # Priority 4: CSV/TSV count matrix (count, dge, expression, ct, etc.)
    csv_patterns = ["*counts*", "*dge*", "*expression*", "*UMI*", "*ct*"]
    csv_files = []
    for pat in csv_patterns:
        csv_files.extend(raw_dir.glob(f"{pat}.csv*"))
        csv_files.extend(raw_dir.glob(f"{pat}.tsv*"))
        csv_files.extend(raw_dir.glob(f"{pat}.txt*"))
    # Exclude summary/metadata files — only keep large count matrices
    csv_files = [f for f in csv_files if "summary" not in f.name.lower() and "readme" not in f.name.lower()]
    # Also exclude RDS files that glob matched
    csv_files = [f for f in csv_files if not f.name.endswith(".rds") and not f.name.endswith(".rds.gz")]
    if csv_files:
        if len(csv_files) == 1:
            return _load_csv_counts(csv_files[0])
        # Multiple count matrices — concatenate them
        logger.info("Found %d count matrices, concatenating...", len(csv_files))
        adatas = []
        for f in sorted(csv_files):
            a = _load_csv_counts(f)
            # Tag cells with source file for traceability
            a.obs["source_file"] = f.name
            adatas.append(a)
        return ad.concat(adatas, join="outer")

    # Priority 5: Any CSV/TSV/TXT file as a last resort (skip metadata-sized files)
    all_tabular = list(raw_dir.glob("*.csv*")) + list(raw_dir.glob("*.tsv*")) + list(raw_dir.glob("*.txt*"))
    all_tabular = [f for f in all_tabular if f.stat().st_size > 50_000]  # skip tiny metadata files
    all_tabular = [f for f in all_tabular if "readme" not in f.name.lower() and "umap" not in f.name.lower()]
    if all_tabular:
        logger.info("No named pattern matched; trying largest tabular file: %s", all_tabular[0].name)
        # Sort by size descending and try the largest
        all_tabular.sort(key=lambda p: p.stat().st_size, reverse=True)
        return _load_csv_counts(all_tabular[0])

    msg = f"No recognized data format found in {raw_dir}"
    raise FileNotFoundError(msg)


def _load_10x_mtx(mtx_dir: Path) -> ad.AnnData:
    """Load 10X-style MTX directory."""
    logger.info("Loading 10X MTX from %s", mtx_dir)

    mtx_path = _find_file(mtx_dir, "matrix.mtx")
    barcodes_path = _find_file(mtx_dir, "barcodes.tsv")
    features_path = _find_file(mtx_dir, "features.tsv", "genes.tsv")

    mat = scipy.io.mmread(str(mtx_path)).T.tocsr()

    barcodes = _read_tsv_column(barcodes_path, col=0)
    features_df = _read_features(features_path)

    adata = ad.AnnData(
        X=mat,
        obs=pd.DataFrame(index=barcodes),
        var=features_df,
    )
    return adata


def _load_csv_counts(path: Path) -> ad.AnnData:
    """Load a CSV/TSV count matrix (genes × cells or cells × genes)."""
    logger.info("Loading count matrix from %s", path)

    # DGE and txt files are typically tab-separated; CSV files use comma
    if ".csv" in path.name:
        sep = ","
    else:
        sep = "\t"
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(fh, sep=sep, index_col=0)

    # Heuristic: if rows >> cols, assume genes × cells → transpose
    if df.shape[0] > df.shape[1] * 2:
        df = df.T

    return ad.AnnData(
        X=scipy.sparse.csr_matrix(df.values),
        obs=pd.DataFrame(index=df.index),
        var=pd.DataFrame(index=df.columns),
    )


def _find_file(directory: Path, *patterns: str) -> Path:
    """Find a file matching any of the patterns (with optional .gz)."""
    for pattern in patterns:
        for suffix in ["", ".gz"]:
            matches = list(directory.glob(f"*{pattern}{suffix}"))
            if matches:
                return matches[0]
    msg = f"Could not find {patterns} in {directory}"
    raise FileNotFoundError(msg)


def _read_tsv_column(path: Path, col: int = 0) -> list[str]:
    """Read a single column from a TSV file (supports .gz)."""
    opener = gzip.open if path.name.endswith(".gz") else open
    values = []
    with opener(path, "rt") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            if parts:
                values.append(parts[col])
    return values


def _read_features(path: Path) -> pd.DataFrame:
    """Read features/genes TSV into a DataFrame."""
    opener = gzip.open if path.name.endswith(".gz") else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            parts = line.strip().split("\t")
            rows.append(parts)

    if rows and len(rows[0]) >= 2:
        gene_ids = [r[0] for r in rows]
        gene_names = [r[1] if len(r) > 1 else r[0] for r in rows]
        df = pd.DataFrame({"gene_ids": gene_ids, "gene_name": gene_names})
        df.index = df["gene_ids"]
    else:
        gene_ids = [r[0] for r in rows]
        df = pd.DataFrame(index=gene_ids)

    return df
