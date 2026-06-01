"""Unit tests for the Q4.5 trigger-element scan (motif_scan + candidates).

Pure-synthetic: consensus PWMs over toy sequences and toy peaks/RepeatMasker
intervals, so no genome FASTA, pysam, or network is needed.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from src.cis_regulatory.candidates import nominate_trigger_elements
from src.cis_regulatory.crossspecies import (
    LiftResult,
    classify_lift,
    convergence_verdict,
    lift_element,
    resolve_call,
)
from src.cis_regulatory.motif_scan import (
    Motif,
    best_relative_score,
    build_background,
    calibrate_thresholds,
    load_motifs,
    scan_motifs,
)

_BASE_IDX = {"A": 0, "C": 1, "G": 2, "T": 3}


def _consensus_motif(name: str, consensus: str) -> Motif:
    """A PWM that strongly prefers ``consensus`` (rel score 1.0 on a match)."""
    counts = np.full((4, len(consensus)), 1, dtype=int)
    for j, base in enumerate(consensus):
        counts[_BASE_IDX[base], j] = 100
    return Motif(matrix_id=name, name=name, counts=counts)


def test_load_real_motifs():
    """The shipped JASPAR config parses into the expected decidual factors."""
    motifs = load_motifs("configs/decidual_motifs.jaspar")
    names = {m.name for m in motifs}
    assert len(motifs) == 8
    assert {"PGR", "NR3C1", "Foxo1", "CEBPB", "GATA2", "Sox17", "ESR1", "Hoxa11"} == names
    # Every motif has a positive length and finite score bounds.
    for m in motifs:
        assert m.length == m.counts.shape[1] > 0
        assert m.max_total > m.min_total


def test_pwm_detects_consensus_and_rejects_mismatch():
    motif = _consensus_motif("M", "ACGTACGT")
    assert best_relative_score("GGGACGTACGTGGG", motif) >= 0.99
    assert best_relative_score("TTTTTTTTTTTT", motif) < 0.85


def test_pwm_detects_reverse_strand():
    motif = _consensus_motif("M", "ACGTACGT")
    # Reverse complement of ACGTACGT is ACGTACGT's revcomp = ACGTACGT -> use
    # a non-palindromic consensus to make the strand test meaningful.
    motif = _consensus_motif("M", "AAAGGGCT")
    revcomp = "AGCCCTTT"  # revcomp of AAAGGGCT
    assert best_relative_score("TT" + revcomp + "TT", motif) >= 0.99


def test_scan_motifs_table_shape_and_hits():
    motifs = [_consensus_motif("PGR", "ACGTACGT"), _consensus_motif("FOXO1", "TTGTTTAC")]
    seqs = {"p1": "GGACGTACGTGG", "p2": "AAAAAAAAAAAA"}
    out = scan_motifs(seqs, motifs)
    assert len(out) == 4  # 2 seqs x 2 motifs
    p1_pgr = out[(out["name"] == "p1") & (out["motif"] == "PGR")].iloc[0]
    assert p1_pgr["hit"]
    p2_pgr = out[(out["name"] == "p2") & (out["motif"] == "PGR")].iloc[0]
    assert not p2_pgr["hit"]


def _peaks() -> pd.DataFrame:
    starts = [1500, 2000, 2500, 50000]
    names = ["pA", "pB", "pC", "pD"]
    return pd.DataFrame(
        {
            "chrom": ["chr1"] * 4,
            "start": starts,
            "end": [s + 60 for s in starts],
            "name": names,
            "score": [40, 30, 20, 10],
            "width": [60] * 4,
        }
    )


def _rmsk() -> pd.DataFrame:
    """MER20 over pA, an Alu over pB, MER20 over the far pD; pC repeat-free."""
    return pd.DataFrame(
        {
            "chrom": ["chr1", "chr1", "chr1"],
            "start": [1450, 1950, 49950],
            "end": [1600, 2100, 50100],
            "repName": ["MER20", "AluY", "MER20"],
            "repClass": ["DNA", "SINE", "DNA"],
            "repFamily": ["hAT-Charlie", "Alu", "hAT-Charlie"],
        }
    )


def _windows() -> pd.DataFrame:
    return pd.DataFrame({"chrom": ["chr1"], "start": [1000], "end": [3000]})


def test_nominate_gates_on_te_proximity_and_pgr():
    pgr = _consensus_motif("PGR", "ACGTACGT")
    foxo = _consensus_motif("FOXO1", "TTGTTTAC")
    motifs = [pgr, foxo]
    seqs = {
        "pA": "GG" + "ACGTACGT" + "GG" + "TTGTTTAC" + "GG",  # PGR + FOXO -> nominated
        "pB": "GGGGGG" + "TTGTTTAC" + "GG",  # FOXO only, no PGR
        "pC": "ACGTACGT" + "TTGTTTAC",  # both, but not TE-derived -> excluded
        "pD": "ACGTACGT" + "TTGTTTAC",  # both + TE, but far from gene -> excluded
    }
    tss = pd.DataFrame({"chrom": ["chr1"], "tss": [2000], "strand": ["+"], "symbol": ["GENEX"]})

    out = nominate_trigger_elements(_peaks(), _rmsk(), _windows(), seqs, motifs, decidual_tss=tss)

    # Only TE-derived near-gene peaks survive: pA and pB (pC no TE, pD far).
    assert set(out["name"]) == {"pA", "pB"}
    nominated = set(out.loc[out["nominated"], "name"])
    assert nominated == {"pA"}

    a = out[out["name"] == "pA"].iloc[0]
    assert a["has_pgr"]
    assert a["n_motifs"] == 2
    assert a["te_flagged"]  # MER20
    assert a["nearest_gene"] == "GENEX"

    b = out[out["name"] == "pB"].iloc[0]
    assert not b["has_pgr"]


# ---------------------------------------------------------------------------
# Per-motif threshold calibration (motif_scan)
# ---------------------------------------------------------------------------


def test_calibration_is_stricter_for_short_motifs():
    """Short PWMs need a higher relative cutoff for the same background FPR."""
    rng = np.random.default_rng(0)
    bg = ["".join(rng.choice(list("ACGT"), size=200)) for _ in range(40)]
    short = _consensus_motif("SHORT", "ACGT")
    long = _consensus_motif("LONG", "ACGTACGTACGTACGT")
    cuts = calibrate_thresholds([short, long], bg, fpr=0.05)
    # A 4-mer hits random DNA far more easily than a 16-mer, so its
    # background-calibrated cutoff must be higher.
    assert cuts["SHORT"] > cuts["LONG"]


def test_build_background_preserves_length_and_composition():
    seqs = {"p1": "AACCGGTT", "p2": "AAAACCCC"}
    bg = build_background(seqs, per_seq=2, seed=1)
    assert len(bg) == 4
    for s in bg:
        assert len(s) == 8
    # Composition of each shuffle matches one of the inputs.
    assert sorted(bg[0]) in (sorted("AACCGGTT"), sorted("AAAACCCC"))


def test_scan_accepts_per_motif_threshold_mapping():
    motif = _consensus_motif("M", "ACGTACGT")
    seqs = {"p1": "GGACGTACGTGG"}
    # An impossible cutoff for M means no hit even on a perfect match.
    out = scan_motifs(seqs, [motif], threshold={"M": 1.01})
    assert not out.iloc[0]["hit"]
    out2 = scan_motifs(seqs, [motif], threshold={"M": 0.5})
    assert out2.iloc[0]["hit"]


# ---------------------------------------------------------------------------
# Cross-species presence/absence + convergence verdict (crossspecies)
# ---------------------------------------------------------------------------


class _FakeLiftOver:
    """Stand-in for pyliftover.LiftOver: maps preset positions, else nothing."""

    def __init__(self, mapped_points, target_chrom="chrT"):
        # mapped_points: set of hg19 positions that successfully lift.
        self._mapped = set(mapped_points)
        self._target = target_chrom

    def convert_coordinate(self, chrom, pos):
        if pos in self._mapped:
            return [(self._target, pos, "+", 0)]
        return []


def test_classify_lift_thresholds():
    assert classify_lift(0.8) == "PRESENT"
    assert classify_lift(0.5) == "PRESENT"
    assert classify_lift(0.2) == "DEGRADED"
    assert classify_lift(0.05) == "ABSENT"


def test_resolve_call_separates_gap_from_loss():
    # Element does not lift but its gene does -> genuine loss (ABSENT).
    loss = LiftResult("mouse", 0.0, "", 0, 0, gene_lifts=True)
    assert resolve_call(loss) == "ABSENT"
    # Neither element nor gene lifts -> uninformative GAP.
    gap = LiftResult("mouse", 0.0, "", 0, 0, gene_lifts=False)
    assert resolve_call(gap) == "GAP"
    # A well-lifting element is PRESENT regardless of the gene flag.
    present = LiftResult("mouse", 0.9, "chrT", 10, 90, gene_lifts=False)
    assert resolve_call(present) == "PRESENT"


def test_lift_element_fraction_and_gene_flag():
    # Lift every sampled point across [100, 200] -> fraction 1.0, PRESENT.
    pts = set(range(100, 201))
    lo = _FakeLiftOver(pts)
    res = lift_element(lo, "mouse", "chr1", 100, 200, gene_chrom="chr1", gene_pos=150)
    assert res.fraction == 1.0
    assert res.gene_lifts
    assert resolve_call(res) == "PRESENT"
    # An element that lifts nowhere, but whose gene lifts -> ABSENT (loss).
    lo2 = _FakeLiftOver({150})  # only the gene TSS lifts
    res2 = lift_element(lo2, "mouse", "chr1", 100, 200, gene_chrom="chr1", gene_pos=150)
    assert res2.fraction == 0.0
    assert res2.gene_lifts
    assert resolve_call(res2) == "ABSENT"


def test_convergence_verdict_labels_trait_tracking_pattern():
    presence = pd.DataFrame(
        {
            "name": ["conv", "cons", "gone", "gapd"],
            "nominated": [True, True, False, False],
            "bat_carollia_class": ["PRESENT", "PRESENT", "ABSENT", "PRESENT"],
            "mouse_class": ["ABSENT", "PRESENT", "ABSENT", "GAP"],
            "ground_squirrel_class": ["DEGRADED", "PRESENT", "ABSENT", "GAP"],
        }
    )
    out = convergence_verdict(presence)
    v = dict(zip(out["name"], out["verdict"], strict=True))
    assert v["conv"] == "CONVERGENT"  # bat present, rodents lost
    assert v["cons"] == "CONSERVED"  # present everywhere
    assert v["gone"] == "ABSENT"  # nowhere
    # bat present but rodents are GAP (uninformative) -> not CONVERGENT.
    assert v["gapd"] == "MIXED"
    # CONVERGENT sorts first.
    assert out.iloc[0]["verdict"] == "CONVERGENT"


def test_crossspecies_presence_table_shape():
    candidates = pd.DataFrame(
        {
            "name": ["e1"],
            "chrom": ["chr1"],
            "start": [100],
            "end": [200],
            "nearest_gene": ["GENEX"],
            "te_name": ["MER41B"],
            "te_family": ["ERV1"],
            "te_flagged": [True],
            "nominated": [True],
        }
    )
    gene_tss = pd.DataFrame({"symbol": ["GENEX"], "chrom": ["chr1"], "tss": [150]})
    import src.cis_regulatory.crossspecies as cs

    fake_lift = _FakeLiftOver(set(range(100, 201)))
    orig = cs._load_liftover
    cs._load_liftover = lambda path: fake_lift  # type: ignore[assignment]
    try:
        out = cs.crossspecies_presence(candidates, {"mouse": "x.chain"}, gene_tss=gene_tss)
    finally:
        cs._load_liftover = orig  # type: ignore[assignment]
    assert list(out["name"]) == ["e1"]
    assert out.iloc[0]["mouse_class"] == "PRESENT"
    assert out.iloc[0]["mouse_frac"] == 1.0
    assert out.iloc[0]["nearest_gene"] == "GENEX"
