"""One-off ingest script for GSE104721 (Sato 2018, EOGT in-vitro decidualization).

Why standalone: the dataset ships a single ``.xls`` with a manually-laid-out
header (Ensembl ID + GeneName + GeneDescription + 12 count columns + 12 TPM
columns) that does not fit the generic ``_load_csv_counts`` path in
``src/ingest/anndata_writer.py``. Rather than complicate the generic loader
for one outlier, this script writes a clean ``results/processed/GSE104721.h5ad``
that downstream QC + scoring can consume unchanged.

Usage:
    python scripts/ingest_gse104721.py

Selects the **siRNA-NT (control)** half of the experiment (6 samples,
3 biological replicates × 2 time-points: NTD0 vs NTD4), discards
HTSeq-stats rows, deduplicates gene symbols, and writes:

    results/processed/GSE104721.h5ad   shape ≈ (6, ~30k)
"""

from __future__ import annotations

import gzip
import logging
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "results/raw/GSE104721/GSE104721_EOGT_counts_TPM.xls.gz"
OUT = ROOT / "results/processed/GSE104721.h5ad"

# siRNA-NT counts columns (control arm: untreated cells in each donor).
# Excludes EOGT-knockdown samples to keep a clean Day-0 → Day-4 contrast.
NT_SAMPLES = [
    "NTD0-sample1",
    "NTD4-sample1",
    "NTD0-sample2",
    "NTD4-sample2",
    "NTD0-sample3",
    "NTD4-sample3",
]


def main() -> None:
    if not SRC.exists():
        msg = (
            f"Missing {SRC}. Download with:\n"
            "  curl -L -o results/raw/GSE104721/GSE104721_EOGT_counts_TPM.xls.gz "
            "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE104nnn/GSE104721/suppl/"
            "GSE104721_EOGT_counts_TPM.xls.gz"
        )
        raise FileNotFoundError(msg)

    tmp = ROOT / "results/raw/GSE104721/GSE104721_EOGT_counts_TPM.xls"
    if not tmp.exists():
        with gzip.open(SRC, "rb") as fh:
            tmp.write_bytes(fh.read())
    logger.info("Reading %s", tmp)

    df = pd.read_excel(tmp, header=1)  # header row 0 is "Counts" / "TPM" labels

    # The first three columns are the gene-level annotation; the next 12 are
    # raw counts (NTD0/NTD4/EOGTD0/EOGTD4 × 3 biological replicates); the last
    # 12 are TPM. We only want NT counts.
    keep_cols = ["Ensembl.GeneID", "GeneName", "GeneDescription", *NT_SAMPLES]
    missing = [c for c in keep_cols if c not in df.columns]
    if missing:
        msg = f"Unexpected column layout — missing {missing}; got {list(df.columns)[:20]}"
        raise ValueError(msg)
    df = df[keep_cols].copy()

    # Drop HTSeq stats rows ("__alignment_not_unique", "__ambiguous", ...).
    df = df[~df["Ensembl.GeneID"].astype(str).str.startswith("__")]
    df = df.dropna(subset=["GeneName"])

    # Collapse duplicate gene symbols by summing counts (rare; e.g. PAR genes).
    counts = df.set_index("GeneName")[NT_SAMPLES].astype(float)
    counts = counts.groupby(level=0).sum()
    logger.info("Counts matrix: %s genes × %s samples", *counts.shape)

    # Build samples × genes AnnData.
    X = scipy.sparse.csr_matrix(counts.values.T.astype(np.float32))
    obs = pd.DataFrame(
        index=NT_SAMPLES,
        data={
            "species": "human",
            "assay": "bulk-RNA-seq",
            "cycle_stage": ["undifferentiated", "decidualized_d4"] * 3,
            "cell_type": "endometrial_stromal",
            "cell_state": "in_vitro_culture",
            "donor": ["donor1", "donor1", "donor2", "donor2", "donor3", "donor3"],
            "sample": NT_SAMPLES,
            "treatment": "siRNA-NT",
            "decidualization_day": [0, 4] * 3,
        },
    )
    var = pd.DataFrame(index=counts.index)
    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.var_names_make_unique()

    OUT.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(OUT)
    logger.info("Wrote %s  shape=%s", OUT, adata.shape)


if __name__ == "__main__":
    main()
