"""Scan peak sequences for decidual transcription-factor motifs (Q4.5).

The Q4.5 "trigger element" idea is that spontaneous decidualization is
conferred by a *specific* TE-derived *cis*-regulatory element carrying the
progesterone-response machinery (a PGR / NR3C1 site) plus the cooperating
decidual factors (FOXO1, CEBPB, HOXA11, GATA2, SOX17, ESR1). To nominate
such elements objectively from the human peak set we need to ask, for each
candidate enhancer, *which* of those factor motifs it carries.

This module is deliberately pip-pure-python: it parses JASPAR-format
position frequency matrices (``configs/decidual_motifs.jaspar``) and scans
sequences with a relative-log-odds PWM scan (the classic min/max-normalised
score used by MATCH / Biopython), so there is no compiled MEME/FIMO
dependency. The only optional binary-backed import is :mod:`pysam`, used
solely for indexed FASTA region extraction and imported lazily.
"""

from __future__ import annotations

import re
from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

# Row order for all matrices in this module.
_BASES = ("A", "C", "G", "T")
_BASE_IDX = {b: i for i, b in enumerate(_BASES)}
_COMPLEMENT = str.maketrans("ACGT", "TGCA")

# One JASPAR record: a header ``>ID  NAME`` followed by four bracketed rows
# (A/C/G/T), each of which may wrap across several physical lines.
_RECORD_RE = re.compile(r">(\S+)\s+(\S+)(.*?)(?=>|\Z)", re.DOTALL)
_ROW_RE = re.compile(r"([ACGT])\s*\[([^\]]+)\]")


@dataclass
class Motif:
    """A position weight matrix with precomputed scan bounds.

    ``counts`` is a ``4 x L`` integer frequency matrix (rows A, C, G, T).
    ``log_odds`` is the natural-log odds vs a uniform background with a
    pseudocount; ``min_total`` / ``max_total`` are the per-position score
    extrema summed across the motif, used to normalise a window score to a
    ``[0, 1]`` relative score independent of motif length.
    """

    matrix_id: str
    name: str
    counts: np.ndarray
    log_odds: np.ndarray = field(init=False)
    min_total: float = field(init=False)
    max_total: float = field(init=False)

    def __post_init__(self) -> None:
        counts = np.asarray(self.counts, dtype=float)
        col_tot = counts.sum(axis=0, keepdims=True)
        # +0.25 pseudocount keeps zero cells finite; background is uniform.
        freq = (counts + 0.25) / (col_tot + 1.0)
        self.log_odds = np.log(freq / 0.25)
        self.min_total = float(self.log_odds.min(axis=0).sum())
        self.max_total = float(self.log_odds.max(axis=0).sum())

    @property
    def length(self) -> int:
        return int(self.counts.shape[1])


def parse_jaspar(text: str) -> list[Motif]:
    """Parse JASPAR-format PFM text (possibly many records) into motifs."""
    motifs: list[Motif] = []
    for matrix_id, name, body in _RECORD_RE.findall(text):
        rows: dict[str, list[int]] = {}
        for base, nums in _ROW_RE.findall(body):
            rows[base] = [int(float(x)) for x in nums.split()]
        if set(rows) != set(_BASES):
            continue
        counts = np.array([rows[b] for b in _BASES], dtype=int)
        motifs.append(Motif(matrix_id=matrix_id, name=name, counts=counts))
    return motifs


def load_motifs(path: str | Path) -> list[Motif]:
    """Load JASPAR PFMs from a file into :class:`Motif` objects."""
    return parse_jaspar(Path(path).read_text())


def _encode(seq: str) -> np.ndarray:
    """Map a sequence to base indices (0-3); non-ACGT becomes -1."""
    arr = np.full(len(seq), -1, dtype=int)
    for base, idx in _BASE_IDX.items():
        arr[np.frombuffer(seq.encode("ascii"), dtype=np.uint8) == ord(base)] = idx
    return arr


def _scan_strand(codes: np.ndarray, motif: Motif) -> float:
    """Best relative score (in ``[0, 1]``) of ``motif`` over one strand."""
    span = motif.length
    if codes.size < span or motif.max_total == motif.min_total:
        return 0.0
    best = 0.0
    rng = motif.max_total - motif.min_total
    for start in range(codes.size - span + 1):
        window = codes[start : start + span]
        if (window < 0).any():  # contains N / non-ACGT — skip
            continue
        score = motif.log_odds[window, np.arange(span)].sum()
        rel = (score - motif.min_total) / rng
        if rel > best:
            best = rel
    return best


def best_relative_score(seq: str, motif: Motif) -> float:
    """Best relative PWM score for ``motif`` on either strand of ``seq``."""
    seq = seq.upper()
    fwd = _scan_strand(_encode(seq), motif)
    rev = _scan_strand(_encode(seq.translate(_COMPLEMENT)[::-1]), motif)
    return max(fwd, rev)


def scan_motifs(
    seqs: dict[str, str],
    motifs: list[Motif],
    threshold: float | Mapping[str, float] = 0.85,
) -> pd.DataFrame:
    """Per-sequence motif hit table.

    For each input sequence and each motif, records the best relative score
    and whether it clears the threshold (a "hit"). ``threshold`` may be a
    single relative cutoff applied to every motif, or a per-motif mapping
    ``name -> cutoff`` (e.g. from :func:`calibrate_thresholds`); motifs
    missing from the mapping fall back to ``0.85``. Returns a long-format
    DataFrame ``name, motif, score, hit`` — one row per (sequence, motif).
    """

    def cut(motif_name: str) -> float:
        if isinstance(threshold, Mapping):
            return float(threshold.get(motif_name, 0.85))
        return float(threshold)

    rows: list[dict[str, object]] = []
    for name, seq in seqs.items():
        for motif in motifs:
            score = best_relative_score(seq, motif)
            rows.append(
                {
                    "name": name,
                    "motif": motif.name,
                    "score": score,
                    "hit": bool(score >= cut(motif.name)),
                }
            )
    return pd.DataFrame(rows)


def shuffle_seq(seq: str, rng: np.random.Generator) -> str:
    """A mononucleotide-composition-preserving shuffle of ``seq``."""
    chars = np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8).copy()
    rng.shuffle(chars)
    return chars.tobytes().decode("ascii")


def calibrate_thresholds(
    motifs: list[Motif],
    background_seqs: list[str],
    fpr: float = 0.01,
) -> dict[str, float]:
    """Per-motif relative-score cutoff at a fixed background false-positive rate.

    Short PWMs clear any single relative cutoff far more often than long
    ones, so a shared threshold is unfair across motifs. For each motif we
    take the best relative score in each composition-matched background
    sequence and set the cutoff to the ``1 - fpr`` quantile of that
    distribution. A "hit" then means the motif scores better than in
    ``(1 - fpr)`` of random sequences of the same base composition — an
    equal false-positive rate for every motif regardless of length.
    """
    cutoffs: dict[str, float] = {}
    q = 1.0 - fpr
    for motif in motifs:
        scores = [best_relative_score(s, motif) for s in background_seqs]
        cutoffs[motif.name] = float(np.quantile(scores, q)) if scores else 0.85
    return cutoffs


def build_background(
    seqs: dict[str, str],
    *,
    per_seq: int = 3,
    seed: int = 0,
) -> list[str]:
    """Composition-matched background by shuffling each input sequence.

    Produces ``per_seq`` shuffles of every sequence in ``seqs``, preserving
    each one's length and base composition — the right null for "is this
    motif really enriched in this element" calibration.
    """
    rng = np.random.default_rng(seed)
    bg: list[str] = []
    for seq in seqs.values():
        for _ in range(per_seq):
            bg.append(shuffle_seq(seq, rng))
    return bg


def extract_peak_seqs(
    peaks: pd.DataFrame,
    fasta_path: str | Path,
) -> dict[str, str]:
    """Extract each peak's reference sequence from an indexed FASTA.

    Uses :mod:`pysam` (lazy import) on a ``.fai``-indexed genome FASTA.
    Chromosome names are tried as-is and with the ``chr`` prefix stripped,
    so a UCSC-style peak table works against an Ensembl-style genome.
    Returns ``name -> uppercase sequence``.
    """
    import pysam

    fa = pysam.FastaFile(str(fasta_path))
    refs = set(fa.references)
    out: dict[str, str] = {}
    for chrom, start, end, name in zip(
        peaks["chrom"], peaks["start"], peaks["end"], peaks["name"], strict=False
    ):
        ref = chrom if chrom in refs else chrom.removeprefix("chr")
        if ref not in refs:
            continue
        out[str(name)] = fa.fetch(ref, int(start), int(end)).upper()
    return out
