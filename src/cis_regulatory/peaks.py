"""Load and summarise GSE61793 MACS peak calls (hg19 BED).

GSE61793 ships BED5 files (``chrom  start  end  name  score``) per assay.
These loaders normalise them to the bioframe column convention
(``chrom``, ``start``, ``end`` plus ``name``/``score``) and provide a
light QC summary, since the deposit is already peak-called — there is no
read-level QC to run, only peak-level sanity stats.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import bioframe as bf
import numpy as np
import pandas as pd

# Series-level supplementary files keyed by a short assay label.
PEAK_FILES: dict[str, str] = {
    "h3k27ac": "GSE61793_h3k27ac_peaks.bed.gz",
    "h3k4me3": "GSE61793_H3K4me3_peaks.bed.gz",
    "dnasei": "GSE61793_DNaseI_80_Union_hg19_.bed.gz",
}

_BED5_COLS = ["chrom", "start", "end", "name", "score"]


def load_peaks(path: str | Path) -> pd.DataFrame:
    """Read a GSE61793 BED5 peak file into a bioframe-style DataFrame.

    The result is sorted, restricted to canonical assembled chromosomes
    (``chr1``..``chr22``, ``chrX``, ``chrY``), and carries an integer
    ``width`` column. ``start``/``end`` keep the BED half-open convention.
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(
            fh,
            sep="\t",
            header=None,
            comment="#",
            usecols=[0, 1, 2, 3, 4],
            names=_BED5_COLS,
            dtype={"chrom": str, "start": int, "end": int, "name": str, "score": float},
        )
    canonical = {f"chr{c}" for c in [*range(1, 23), "X", "Y"]}
    df = df[df["chrom"].isin(canonical)].copy()
    df = bf.sort_bedframe(df)
    df["width"] = df["end"] - df["start"]
    return df.reset_index(drop=True)


def peak_qc(peaks: pd.DataFrame, assay: str) -> dict[str, object]:
    """Summary statistics for one assay's peak set.

    Returns a flat dict (n_peaks, total/median/mean width, chromosome
    count) suitable for tabulating across assays.
    """
    widths = peaks["width"].to_numpy()
    return {
        "assay": assay,
        "n_peaks": int(len(peaks)),
        "n_chroms": int(peaks["chrom"].nunique()),
        "total_bp": int(widths.sum()),
        "median_width": float(np.median(widths)) if len(widths) else float("nan"),
        "mean_width": float(widths.mean()) if len(widths) else float("nan"),
    }


def load_all_peaks(raw_dir: str | Path) -> dict[str, pd.DataFrame]:
    """Load every available GSE61793 peak assay found under ``raw_dir``.

    Missing files are skipped (not all assays are required for every
    analysis), so the caller can run with whatever is on disk.
    """
    raw_dir = Path(raw_dir)
    out: dict[str, pd.DataFrame] = {}
    for assay, fname in PEAK_FILES.items():
        fpath = raw_dir / fname
        if fpath.exists():
            out[assay] = load_peaks(fpath)
    return out
