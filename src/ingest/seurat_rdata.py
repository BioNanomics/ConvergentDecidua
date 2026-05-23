"""Convert Seurat ``.RData`` files to AnnData via an Rscript bridge.

Requires ``R`` + ``Seurat`` + ``Matrix`` on PATH (see ``docs/REPRODUCE.md`` and
``docker/Dockerfile``). The R script writes a 10X-style MTX directory plus an
``obs.csv``; this module loads them into AnnData and (optionally) concatenates
across multiple RData files in a directory.
"""

from __future__ import annotations

import logging
import shutil
import subprocess
from pathlib import Path

import anndata as ad
import pandas as pd
import scipy.io
import scipy.sparse

logger = logging.getLogger(__name__)

R_SCRIPT = Path(__file__).resolve().parents[2] / "scripts" / "rdata_to_h5ad.R"


def rdata_to_anndata(
    rdata_path: Path,
    work_dir: Path,
    *,
    object_name: str | None = None,
    rscript: str = "Rscript",
) -> ad.AnnData:
    """Convert a single Seurat ``.RData(.gz)`` file to AnnData."""
    if shutil.which(rscript) is None:
        msg = (
            f"`{rscript}` not found on PATH. Install R + Seurat (see "
            "docs/REPRODUCE.md), or run this stage inside the Docker image."
        )
        raise RuntimeError(msg)

    work_dir.mkdir(parents=True, exist_ok=True)
    cmd = [rscript, str(R_SCRIPT), str(rdata_path), str(work_dir)]
    if object_name:
        cmd.append(object_name)
    logger.info("Running R bridge: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return _load_mtx_with_obs(work_dir)


def rdata_dir_to_anndata(
    raw_dir: Path,
    work_root: Path,
    *,
    include: list[str] | None = None,
    object_name: str | None = None,
    rscript: str = "Rscript",
) -> ad.AnnData:
    """Convert every ``.RData(.gz)`` in ``raw_dir`` and concatenate.

    Parameters
    ----------
    raw_dir
        Directory containing one or more ``.RData(.gz)`` files.
    work_root
        Directory under which a sub-folder is created per input file.
    include
        Optional list of filename substrings; only matching files are
        converted. Useful for selecting (e.g.) ``["UE_DSC"]`` from a
        multi-tissue dump.
    """
    files = sorted(
        p for p in raw_dir.iterdir() if p.name.endswith((".RData", ".RData.gz"))
    )
    if include:
        files = [p for p in files if any(tag in p.name for tag in include)]
    if not files:
        msg = f"No matching .RData files in {raw_dir} (include={include})"
        raise FileNotFoundError(msg)

    adatas: list[ad.AnnData] = []
    for f in files:
        sub = work_root / f.name.replace(".RData.gz", "").replace(".RData", "")
        logger.info("Converting %s -> %s", f.name, sub)
        a = rdata_to_anndata(f, sub, object_name=object_name, rscript=rscript)
        a.obs["source_file"] = f.name
        adatas.append(a)

    if len(adatas) == 1:
        return adatas[0]
    return ad.concat(adatas, join="outer", merge="unique")


def _load_mtx_with_obs(mtx_dir: Path) -> ad.AnnData:
    mat = scipy.io.mmread(str(mtx_dir / "matrix.mtx")).T.tocsr()
    barcodes = [
        ln.strip()
        for ln in (mtx_dir / "barcodes.tsv").read_text().splitlines()
        if ln.strip()
    ]
    features = pd.read_csv(
        mtx_dir / "features.tsv",
        sep="\t",
        header=None,
        names=["gene_id", "gene_symbol"],
    )
    obs = pd.read_csv(mtx_dir / "obs.csv")
    obs.index = obs["barcode"].astype(str)
    obs = obs.reindex(barcodes)
    var = pd.DataFrame(
        {"gene_symbol": features["gene_symbol"].values},
        index=features["gene_id"].astype(str).values,
    )
    return ad.AnnData(X=scipy.sparse.csr_matrix(mat), obs=obs, var=var)
