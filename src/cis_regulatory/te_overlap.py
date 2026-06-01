"""Overlap GSE61793 peaks with RepeatMasker TEs (Lynch/Wagner test).

The Lynch/Wagner model holds that the decidual *cis*-regulatory landscape
was substantially built from ancient transposable elements — notably the
MER20 (hAT-Charlie, DNA) and MER41 (ERV1, LTR) families. This module
quantifies, for each assay's peak set, the fraction of peaks that overlap
a TE at all, the fraction by TE class, and the fraction overlapping the
flagged MER20 / MER41 families specifically.

Coordinates are UCSC hg19 / BED half-open throughout (both GSE61793 and
RepeatMasker), so no liftover is needed. Overlap is computed with
:mod:`bioframe`.
"""

from __future__ import annotations

import gzip
from pathlib import Path

import bioframe as bf
import pandas as pd
from scipy import stats

# UCSC ``rmsk.txt`` column positions (no header in the file).
_RMSK_USECOLS = [5, 6, 7, 10, 11, 12]
_RMSK_NAMES = ["chrom", "start", "end", "repName", "repClass", "repFamily"]

# RepeatMasker classes that are genuine transposable elements (excludes
# Simple_repeat, Low_complexity, Satellite, rRNA/tRNA/snRNA, Unknown). In
# UCSC hg19 the SVA retroposons are filed under class ``Other`` (not
# ``Retroposon``), so ``Other`` is included to capture them.
TE_CLASSES = ("DNA", "LINE", "SINE", "LTR", "Other", "RC")

# Families Lynch et al. implicated in decidual cis-rewiring.
FLAGGED_FAMILIES = ("MER20", "MER41")


def load_rmsk(path: str | Path) -> pd.DataFrame:
    """Read UCSC ``rmsk.txt.gz`` into a bioframe-style TE table.

    Returns ``chrom, start, end, repName, repClass, repFamily`` restricted
    to canonical assembled chromosomes and sorted. ``repClass`` values
    sometimes carry a ``/`` qualifier (e.g. ``DNA?``); the trailing ``?``
    (RepeatMasker's low-confidence marker) is stripped so class filters
    match.
    """
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(
            fh,
            sep="\t",
            header=None,
            usecols=_RMSK_USECOLS,
            names=_RMSK_NAMES,
            dtype={c: str for c in _RMSK_NAMES} | {"start": int, "end": int},
        )
    canonical = {f"chr{c}" for c in [*range(1, 23), "X", "Y"]}
    df = df[df["chrom"].isin(canonical)].copy()
    df["repClass"] = df["repClass"].str.rstrip("?")
    return bf.sort_bedframe(df).reset_index(drop=True)


def _overlap_stats(peaks: pd.DataFrame, intervals: pd.DataFrame) -> tuple[int, float]:
    """How many of ``peaks`` overlap >=1 row of ``intervals`` (count, frac)."""
    n = len(peaks)
    if n == 0 or len(intervals) == 0:
        return 0, 0.0
    cov = bf.count_overlaps(peaks, intervals[["chrom", "start", "end"]])
    n_hit = int((cov["count"] > 0).sum())
    return n_hit, n_hit / n


def te_enrichment(
    peaks: pd.DataFrame,
    rmsk: pd.DataFrame,
    assay: str,
    flagged_families: tuple[str, ...] = FLAGGED_FAMILIES,
) -> pd.DataFrame:
    """TE-derived fractions for one assay's peak set.

    Returns one row per category with ``n_peaks`` (total in the assay),
    ``n_overlap`` and ``frac_overlap``:

    - ``any_repeat`` — overlaps any RepeatMasker element,
    - ``any_TE`` — overlaps a genuine TE class (:data:`TE_CLASSES`),
    - ``class:<C>`` — overlaps that TE class,
    - ``family:<F>*`` — overlaps a flagged family (``repName`` prefix),
      the direct MER20 / MER41 Lynch-model test.
    """
    te = rmsk[rmsk["repClass"].isin(TE_CLASSES)]
    n_peaks = len(peaks)

    rows: list[dict[str, object]] = []

    def add(category: str, intervals: pd.DataFrame) -> None:
        n_hit, frac = _overlap_stats(peaks, intervals)
        rows.append(
            {
                "assay": assay,
                "category": category,
                "n_peaks": n_peaks,
                "n_overlap": n_hit,
                "frac_overlap": frac,
            }
        )

    add("any_repeat", rmsk)
    add("any_TE", te)
    for cls in TE_CLASSES:
        add(f"class:{cls}", te[te["repClass"] == cls])
    for fam in flagged_families:
        add(f"family:{fam}*", rmsk[rmsk["repName"].str.startswith(fam)])

    return pd.DataFrame(rows)


def te_enrichment_all(
    peaks_by_assay: dict[str, pd.DataFrame],
    rmsk: pd.DataFrame,
    flagged_families: tuple[str, ...] = FLAGGED_FAMILIES,
) -> pd.DataFrame:
    """Run :func:`te_enrichment` across every assay and concatenate."""
    frames = [
        te_enrichment(peaks, rmsk, assay, flagged_families)
        for assay, peaks in peaks_by_assay.items()
    ]
    return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()


def _flag_intervals(rmsk: pd.DataFrame, family: str) -> pd.DataFrame:
    """RepeatMasker rows whose ``repName`` starts with ``family``."""
    return rmsk[rmsk["repName"].str.startswith(family)][["chrom", "start", "end"]]


def near_gene_te_enrichment(
    peaks: pd.DataFrame,
    rmsk: pd.DataFrame,
    gene_windows: pd.DataFrame,
    assay: str,
    flagged_families: tuple[str, ...] = FLAGGED_FAMILIES,
) -> pd.DataFrame:
    """Lynch test: is a flagged TE family enriched in decidual-gene peaks?

    Partitions ``peaks`` into those overlapping ``gene_windows`` (decidual
    regulatory neighborhoods) vs the rest, then for each flagged family
    asks whether family overlap is **more frequent among near-gene
    peaks** via a one-sided Fisher's exact test (alternative="greater").
    A significant enrichment of MER20 / MER41 in decidual-gene enhancers
    is the direct prediction of the Lynch/Wagner cis-rewiring model.

    Returns one row per family with the 2x2 counts, the near/far overlap
    fractions, the odds ratio and the one-sided p-value.
    """
    n_peaks = len(peaks)
    if n_peaks == 0 or len(gene_windows) == 0:
        return pd.DataFrame()

    near_mask = bf.count_overlaps(peaks, gene_windows[["chrom", "start", "end"]])["count"] > 0
    near = peaks[near_mask.to_numpy()]
    far = peaks[~near_mask.to_numpy()]

    rows: list[dict[str, object]] = []
    for fam in flagged_families:
        fam_iv = _flag_intervals(rmsk, fam)
        near_hit, near_frac = _overlap_stats(near, fam_iv)
        far_hit, far_frac = _overlap_stats(far, fam_iv)
        table = [[near_hit, len(near) - near_hit], [far_hit, len(far) - far_hit]]
        odds, pval = stats.fisher_exact(table, alternative="greater")
        rows.append(
            {
                "assay": assay,
                "family": f"{fam}*",
                "n_near": int(len(near)),
                "n_far": int(len(far)),
                "near_hit": int(near_hit),
                "far_hit": int(far_hit),
                "near_frac": near_frac,
                "far_frac": far_frac,
                "odds_ratio": float(odds),
                "fisher_p_greater": float(pval),
            }
        )
    return pd.DataFrame(rows)
