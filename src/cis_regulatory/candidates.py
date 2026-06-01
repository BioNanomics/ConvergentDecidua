"""Nominate candidate TE-derived decidual "trigger" elements (Q4.5).

The convergence hypothesis needs a *specific* candidate locus to test
across species, and the Q3.4 expression-correlation ranking cannot supply
one: the canonical initiators (PGR, HAND2, CEBPB, HOXA11, GATA2, SOX17)
were dropped at highly-variable-gene selection and the ranking is
structurally blind to lowly-expressed, ligand-gated master switches. So we
nominate candidates from *sequence* instead.

A candidate trigger element is a human H3K27ac enhancer peak that is
simultaneously:

1. **TE-derived** — overlaps a genuine RepeatMasker TE class (the
   Lynch/Wagner substrate; MER20 / MER41 families flagged),
2. **decidual-proximal** — within the ±window decidual-gene neighborhoods,
3. **progesterone-responsive** — carries a PGR / NR3C1 motif plus at least
   one more cooperating decidual factor motif (FOXO1, CEBPB, HOXA11,
   GATA2, SOX17, ESR1).

This is the unbiased "scan" half of the scan-then-drill design: it ranks
elements by how trigger-like they are so the expensive cross-species
drill-down runs only on the top hits.
"""

from __future__ import annotations

from collections.abc import Mapping

import bioframe as bf
import pandas as pd

from src.cis_regulatory.motif_scan import Motif, scan_motifs
from src.cis_regulatory.te_overlap import FLAGGED_FAMILIES, TE_CLASSES

# Motifs that satisfy the progesterone-response requirement (gate condition).
PGR_MOTIFS = ("PGR", "NR3C1")

_TE_COLS = ["chrom", "start", "end", "repName", "repClass", "repFamily"]


def _annotate_te(peaks: pd.DataFrame, rmsk: pd.DataFrame) -> pd.DataFrame:
    """Annotate each peak with its overlapping TE, preferring flagged families.

    Inner overlap against genuine TE classes, so the result keeps only
    TE-derived peaks. When a peak overlaps several TEs the flagged
    (MER20/MER41) one wins, else the first. Adds ``te_name``, ``te_class``,
    ``te_family`` and a boolean ``te_flagged``.
    """
    te = rmsk[rmsk["repClass"].isin(TE_CLASSES)][_TE_COLS]
    ov = bf.overlap(peaks, te, how="inner", suffixes=("", "_te"))
    if ov.empty:
        return ov
    ov["te_flagged"] = ov["repName_te"].apply(lambda s: str(s).startswith(FLAGGED_FAMILIES))
    ov = ov.sort_values("te_flagged", ascending=False)
    best = ov.drop_duplicates("name", keep="first").copy()
    return best.rename(
        columns={"repName_te": "te_name", "repClass_te": "te_class", "repFamily_te": "te_family"}
    )


def _nearest_gene(peaks: pd.DataFrame, decidual_tss: pd.DataFrame) -> pd.DataFrame:
    """Nearest decidual-gene symbol and distance for each peak (bioframe closest)."""
    genes = decidual_tss.assign(start=decidual_tss["tss"], end=decidual_tss["tss"] + 1)
    cl = bf.closest(peaks, genes[["chrom", "start", "end", "symbol"]], suffixes=("", "_g"))
    return pd.DataFrame(
        {
            "name": cl["name"].to_numpy(),
            "nearest_gene": cl["symbol_g"].to_numpy(),
            "gene_dist": cl["distance"].to_numpy(),
        }
    )


def nominate_trigger_elements(
    peaks: pd.DataFrame,
    rmsk: pd.DataFrame,
    gene_windows: pd.DataFrame,
    seqs: dict[str, str],
    motifs: list[Motif],
    *,
    decidual_tss: pd.DataFrame | None = None,
    min_motifs: int = 2,
    threshold: float | Mapping[str, float] = 0.85,
    pgr_motifs: tuple[str, ...] = PGR_MOTIFS,
) -> pd.DataFrame:
    """Rank TE-derived, decidual-proximal enhancers by trigger-likeness.

    ``seqs`` maps each peak's ``name`` to its reference sequence (from
    :func:`src.cis_regulatory.motif_scan.extract_peak_seqs`). Returns one
    row per TE-derived near-gene peak with its TE identity, the decidual
    motifs it carries, and a boolean ``nominated`` (passes the gate:
    ``>= min_motifs`` distinct motifs including a PGR/NR3C1 site). Rows are
    sorted best-first; an empty frame means no near-gene TE-derived peaks.
    """
    peaks = peaks.reset_index(drop=True).copy()
    if peaks["name"].duplicated().any():
        peaks["name"] = [f"peak_{i}" for i in range(len(peaks))]

    near_mask = bf.count_overlaps(peaks, gene_windows[["chrom", "start", "end"]])["count"] > 0
    peaks = peaks[near_mask.to_numpy()].copy()
    if peaks.empty:
        return pd.DataFrame()

    ann = _annotate_te(peaks, rmsk)
    if ann.empty:
        return pd.DataFrame()

    sub_seqs = {n: seqs[n] for n in ann["name"] if n in seqs}
    hits = scan_motifs(sub_seqs, motifs, threshold=threshold)
    hits = hits[hits["hit"]]
    by_peak = hits.groupby("name")["motif"].agg(lambda s: sorted(set(s)))

    ann["motifs"] = ann["name"].map(by_peak).apply(lambda v: v if isinstance(v, list) else [])
    ann["n_motifs"] = ann["motifs"].apply(len)
    pgr_set = set(pgr_motifs)
    ann["has_pgr"] = ann["motifs"].apply(lambda v: bool(pgr_set.intersection(v)))
    ann["motifs"] = ann["motifs"].apply(lambda v: ",".join(v))
    ann["nominated"] = (ann["n_motifs"] >= min_motifs) & ann["has_pgr"]

    if decidual_tss is not None and not decidual_tss.empty:
        ann = ann.merge(_nearest_gene(ann, decidual_tss), on="name", how="left")

    sort_cols = ["nominated", "te_flagged", "n_motifs", "score"]
    if "score" not in ann.columns:
        sort_cols.remove("score")
    ann = ann.sort_values(sort_cols, ascending=False).reset_index(drop=True)

    cols = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "width",
        "te_name",
        "te_class",
        "te_family",
        "te_flagged",
        "nearest_gene",
        "gene_dist",
        "n_motifs",
        "has_pgr",
        "motifs",
        "nominated",
    ]
    keep = [c for c in cols if c in ann.columns]
    return ann[keep]
