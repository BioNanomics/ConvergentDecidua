"""Build decidual-gene regulatory neighborhoods from UCSC refGene (hg19).

The Lynch/Wagner TE-rewiring hypothesis is specifically about the
*cis*-regulatory neighborhoods of decidual / progesterone-response genes,
not the genome at large. This module turns the decidual marker gene set
(human symbols from ``configs/markers.yaml``) into hg19 windows around
each transcription start site, so peaks can be partitioned into
"near a decidual gene" vs "elsewhere".
"""

from __future__ import annotations

import gzip
from pathlib import Path

import bioframe as bf
import pandas as pd

# UCSC ``refGene.txt`` column positions (no header).
_REFGENE_USECOLS = [2, 3, 4, 5, 12]
_REFGENE_NAMES = ["chrom", "strand", "txStart", "txEnd", "symbol"]


def load_tss(path: str | Path) -> pd.DataFrame:
    """Read UCSC ``refGene.txt.gz`` into a per-transcript TSS table.

    Returns ``chrom, tss, strand, symbol`` for canonical chromosomes. The
    TSS is ``txStart`` on the ``+`` strand and ``txEnd`` on the ``-``
    strand. Multiple transcripts per gene are kept (collapsed later when
    windows are merged).
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(
            fh,
            sep="\t",
            header=None,
            usecols=_REFGENE_USECOLS,
            names=_REFGENE_NAMES,
            dtype={"chrom": str, "strand": str, "txStart": int, "txEnd": int, "symbol": str},
        )
    canonical = {f"chr{c}" for c in [*range(1, 23), "X", "Y"]}
    df = df[df["chrom"].isin(canonical)].copy()
    df["tss"] = df["txEnd"].where(df["strand"] == "-", df["txStart"])
    return df[["chrom", "tss", "strand", "symbol"]].reset_index(drop=True)


def gene_windows(
    tss: pd.DataFrame,
    symbols: list[str],
    window: int = 50_000,
) -> pd.DataFrame:
    """``±window`` bp hg19 intervals around the TSS of the given symbols.

    Symbol matching is case-insensitive. Overlapping windows (e.g. from
    multiple transcripts of the same gene, or neighbouring genes) are
    merged so each base is counted once. Returns a sorted, merged
    bioframe (``chrom, start, end``); missing symbols are silently
    skipped (reported by the caller via the returned ``symbols`` set).
    """
    wanted = {s.upper() for s in symbols}
    sub = tss[tss["symbol"].str.upper().isin(wanted)].copy()
    if sub.empty:
        return pd.DataFrame(columns=["chrom", "start", "end"])
    sub["start"] = (sub["tss"] - window).clip(lower=0)
    sub["end"] = sub["tss"] + window
    merged = bf.merge(bf.sort_bedframe(sub[["chrom", "start", "end"]]))
    return merged[["chrom", "start", "end"]].reset_index(drop=True)


def matched_symbols(tss: pd.DataFrame, symbols: list[str]) -> set[str]:
    """Subset of ``symbols`` that exist in the refGene table (upper-cased)."""
    present = set(tss["symbol"].str.upper())
    return {s.upper() for s in symbols} & present
