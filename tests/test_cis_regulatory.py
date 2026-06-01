"""Unit tests for src/cis_regulatory (Q4.3) — synthetic peaks + RepeatMasker."""

from __future__ import annotations

import pandas as pd

from src.cis_regulatory.genes import gene_windows, matched_symbols
from src.cis_regulatory.peaks import peak_qc
from src.cis_regulatory.te_overlap import (
    near_gene_te_enrichment,
    te_enrichment,
    te_enrichment_all,
)


def _peaks() -> pd.DataFrame:
    """Four peaks on chr1; widths 100. Peaks at 1000, 5000, 9000, 50000."""
    starts = [1000, 5000, 9000, 50000]
    return pd.DataFrame(
        {
            "chrom": ["chr1"] * 4,
            "start": starts,
            "end": [s + 100 for s in starts],
            "name": [f"p{i}" for i in range(4)],
            "score": [10, 20, 30, 40],
            "width": [100] * 4,
        }
    )


def _rmsk() -> pd.DataFrame:
    """Synthetic RepeatMasker: a SINE over peak0, a MER20 DNA over peak1,
    a MER41 LTR over peak3; peak2 is repeat-free."""
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1"],
            "start": [950, 4950, 49950],
            "end": [1200, 5200, 50200],
            "repName": ["AluY", "MER20", "MER41B"],
            "repClass": ["SINE", "DNA", "LTR"],
            "repFamily": ["Alu", "hAT-Charlie", "ERV1"],
        }
    )


def _tss() -> pd.DataFrame:
    """One decidual gene GENEX with TSS at chr1:5000 (near peak1)."""
    return pd.DataFrame({"chrom": ["chr1"], "tss": [5000], "strand": ["+"], "symbol": ["GENEX"]})


def test_peak_qc_counts():
    qc = peak_qc(_peaks(), "h3k27ac")
    assert qc["assay"] == "h3k27ac"
    assert qc["n_peaks"] == 4
    assert qc["n_chroms"] == 1
    assert qc["median_width"] == 100


def test_te_enrichment_categories_and_fractions():
    out = te_enrichment(_peaks(), _rmsk(), "h3k27ac")
    rows = {r["category"]: r for _, r in out.iterrows()}
    # 3 of 4 peaks overlap a repeat / TE.
    assert rows["any_repeat"]["n_overlap"] == 3
    assert rows["any_TE"]["n_overlap"] == 3
    assert rows["class:SINE"]["n_overlap"] == 1
    assert rows["class:DNA"]["n_overlap"] == 1
    assert rows["class:LTR"]["n_overlap"] == 1
    assert rows["family:MER20*"]["n_overlap"] == 1
    assert rows["family:MER41*"]["n_overlap"] == 1
    assert abs(rows["any_TE"]["frac_overlap"] - 0.75) < 1e-9


def test_te_enrichment_all_concats_assays():
    peaks_by_assay = {"h3k27ac": _peaks(), "h3k4me3": _peaks()}
    out = te_enrichment_all(peaks_by_assay, _rmsk())
    assert set(out["assay"]) == {"h3k27ac", "h3k4me3"}


def test_gene_windows_and_matched_symbols():
    tss = _tss()
    assert matched_symbols(tss, ["GENEX", "ABSENT"]) == {"GENEX"}
    win = gene_windows(tss, ["GENEX"], window=2000)
    assert len(win) == 1
    assert win.iloc[0]["start"] == 3000
    assert win.iloc[0]["end"] == 7000


def test_near_gene_te_enrichment_partitions_peaks():
    # ±2kb window around TSS 5000 → only peak1 (start 5000) is "near".
    win = gene_windows(_tss(), ["GENEX"], window=2000)
    out = near_gene_te_enrichment(_peaks(), _rmsk(), win, "h3k27ac")
    mer20 = out[out["family"] == "MER20*"].iloc[0]
    assert mer20["n_near"] == 1
    assert mer20["near_hit"] == 1  # peak1 carries MER20
    assert mer20["n_far"] == 3
    assert mer20["far_hit"] == 0
    # MER20 perfectly concentrated in the near set → OR is large.
    assert mer20["odds_ratio"] > 1.0
