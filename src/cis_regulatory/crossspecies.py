"""Cross-species presence/absence of candidate trigger elements (Q4.5 Phase B).

The convergence hypothesis predicts a *specific* phylogenetic pattern for a
TE-derived decidual trigger element: **present** in the spontaneously-
decidualizing lineages (human, *Carollia* bat) and **absent or degraded**
in the trait-negative rodents (mouse, ground squirrel). This module is the
"drill" half of the scan-then-drill design.

**Why synteny / liftOver, not raw whole-genome alignment.** An early version
aligned each human element against the target genomes with ``mappy``
(minimap2, ``asm20``). That is *seed-length blind*: across the ~80 My
separating these lineages a ~1.7 kb regulatory element retains no minimizer-
length anchor, so minimap2 returns nothing — even next to deeply conserved
genes (LAMA4, CDKN1A) that are unquestionably present in mouse. Every
element scored ABSENT, giving the test zero discriminating power. UCSC
**liftOver chains** are built from sensitive lastz whole-genome alignments
and are the standard way to ask "where is this human region in mouse"; they
resolve conserved elements that minimap2's ``asm`` presets cannot seed.

For each element we lift its hg19 span through the chain and classify by the
**fraction of the element that lifts** (a conserved element lifts as a
contiguous block; an eroded or absent one lifts partially or not at all). To
separate genuine *element loss* from a mere *assembly/synteny gap*, we also
lift the flanking decidual-gene TSS: if the gene lifts but the element does
not, that is candidate lineage-specific loss; if neither lifts, the region
is simply not resolvable in that assembly and the call is uninformative.

*Carollia* has no UCSC chain (the Bat1K mCarPer1.2 assembly is new), so for
that one trait-positive lineage we fall back to a *sensitive* ``mappy``
search (short-seed preset) anchored on the flanking conserved sequence.

Classification is intentionally coarse (PRESENT / DEGRADED / ABSENT /
GAP) because it feeds a convergence *verdict*, not a base-level claim.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd

# Fraction-lifted thresholds for the presence call.
PRESENT_FRACTION = 0.50
DEGRADED_FRACTION = 0.10

# Number of points sampled across an element span when lifting.
LIFT_POINTS = 20

# Sensitive mappy preset for the chain-less Carollia fallback. "sr" uses
# short-read seeding (small k/w) so a diverged but real core still seeds.
CAROLLIA_PRESET = "sr"
CAROLLIA_PRESENT_COV = 0.30
CAROLLIA_PRESENT_ID = 0.65

# Trait phylogeny for the convergence test. Human is the query/reference
# (always "present" by construction). The other spontaneously-decidualizing
# lineage is the Carollia bat; the trait-negative outgroups are the rodents.
TRAIT_POSITIVE = ("bat_carollia",)
TRAIT_NEGATIVE = ("mouse", "ground_squirrel")


@dataclass
class LiftResult:
    """Outcome of lifting one element span through a chain to one genome."""

    species: str
    fraction: float
    target_chrom: str
    target_start: int
    target_end: int
    gene_lifts: bool

    @property
    def locus(self) -> str:
        if not self.target_chrom:
            return ""
        return f"{self.target_chrom}:{self.target_start}-{self.target_end}"


def _lift_span(lo, chrom: str, start: int, end: int, n: int = LIFT_POINTS):
    """Lift ``n`` evenly-spaced points across a span; summarise the result.

    Returns ``(fraction_lifted, target_chrom, target_start, target_end)``.
    The fraction is the share of sampled points that lift at all; the target
    span is reported only when the lifted points agree on a single target
    chromosome (a coherent syntenic block rather than scattered fragments).
    """
    if end <= start:
        end = start + 1
    pts = [int(start + (end - start) * i / (n - 1)) for i in range(n)]
    mapped = [hit[0] for p in pts if (hit := lo.convert_coordinate(chrom, p))]
    if not mapped:
        return 0.0, "", 0, 0
    fraction = len(mapped) / n
    chroms = {m[0] for m in mapped}
    if len(chroms) == 1:
        positions = [m[1] for m in mapped]
        return fraction, mapped[0][0], min(positions), max(positions)
    return fraction, "", 0, 0


def classify_lift(fraction: float) -> str:
    """PRESENT / DEGRADED / ABSENT from the fraction of an element that lifts."""
    if fraction >= PRESENT_FRACTION:
        return "PRESENT"
    if fraction >= DEGRADED_FRACTION:
        return "DEGRADED"
    return "ABSENT"


def lift_element(
    lo,
    species: str,
    chrom: str,
    start: int,
    end: int,
    *,
    gene_chrom: str | None = None,
    gene_pos: int | None = None,
) -> LiftResult:
    """Lift one element (and, optionally, its flanking gene TSS) to a genome."""
    fraction, t_chrom, t_start, t_end = _lift_span(lo, chrom, start, end)
    gene_lifts = False
    if gene_chrom is not None and gene_pos is not None:
        gene_lifts = bool(lo.convert_coordinate(gene_chrom, int(gene_pos)))
    return LiftResult(
        species=species,
        fraction=fraction,
        target_chrom=t_chrom,
        target_start=t_start,
        target_end=t_end,
        gene_lifts=gene_lifts,
    )


def resolve_call(lift: LiftResult) -> str:
    """Final per-species call, separating loss from an uninformative gap.

    An element that does not lift is only evidence of *loss* when the region
    is otherwise resolvable in that assembly (the flanking gene lifts). If
    neither the element nor its gene lifts, the locus is an assembly/synteny
    **GAP** and the call is uninformative rather than ABSENT.
    """
    base = classify_lift(lift.fraction)
    if base == "ABSENT" and not lift.gene_lifts:
        return "GAP"
    return base


def _carollia_call(seq: str, aligner) -> tuple[str, float, float, str]:
    """Sensitive mappy presence call for the chain-less Carollia genome."""
    if not seq:
        return "GAP", 0.0, 0.0, ""
    qlen = len(seq)
    best = None
    for hit in aligner.map(seq):
        coverage = (hit.q_en - hit.q_st) / qlen if qlen else 0.0
        identity = hit.mlen / hit.blen if hit.blen else 0.0
        score = coverage * identity
        if best is None or score > best[0]:
            best = (score, coverage, identity, f"{hit.ctg}:{hit.r_st}-{hit.r_en}")
    if best is None:
        return "ABSENT", 0.0, 0.0, ""
    _, coverage, identity, locus = best
    if coverage >= CAROLLIA_PRESENT_COV and identity >= CAROLLIA_PRESENT_ID:
        call = "PRESENT"
    elif coverage >= CAROLLIA_PRESENT_COV / 2:
        call = "DEGRADED"
    else:
        call = "ABSENT"
    return call, round(identity, 3), round(coverage, 3), locus


def _load_liftover(chain_path: str | Path):
    """Build a :class:`pyliftover.LiftOver` from a UCSC over.chain(.gz)."""
    from pyliftover import LiftOver

    return LiftOver(str(chain_path))


def _load_carollia_aligner(genome_path: str | Path, preset: str):
    """Build (or load cached ``.mmi``) a sensitive mappy aligner for Carollia."""
    import mappy as mp

    genome_path = Path(genome_path)
    mmi = genome_path.with_suffix(genome_path.suffix + ".mmi")
    if mmi.exists():
        aligner = mp.Aligner(str(mmi), preset=preset)
    else:
        aligner = mp.Aligner(str(genome_path), preset=preset, fn_idx_out=str(mmi))
    if not aligner:
        raise RuntimeError(f"mappy failed to build/load an index for {genome_path}")
    return aligner


def crossspecies_presence(
    candidates: pd.DataFrame,
    chains: dict[str, str | Path],
    *,
    carollia_seqs: dict[str, str] | None = None,
    carollia_genome: str | Path | None = None,
    carollia_preset: str = CAROLLIA_PRESET,
    gene_tss: pd.DataFrame | None = None,
) -> pd.DataFrame:
    """Presence/absence table for each candidate element across target genomes.

    ``candidates`` are human elements (``chrom/start/end/name`` plus carried
    annotation), ``chains`` maps a rodent species label to its hg19->target
    UCSC chain. ``gene_tss`` (``symbol``/``chrom``/``tss``) lets each element
    also lift its nearest decidual gene to separate loss from assembly gaps.
    Carollia is handled by sensitive ``mappy`` when ``carollia_seqs`` and
    ``carollia_genome`` are given. Returns one row per element with, per
    species, ``<sp>_class`` / ``<sp>_frac`` (or ``_ident``/``_cov`` for
    Carollia) / ``<sp>_locus`` — ready for :func:`convergence_verdict`.
    """
    lifters = {sp: _load_liftover(path) for sp, path in chains.items()}
    tss_lookup: dict[str, tuple[str, int]] = {}
    if gene_tss is not None and not gene_tss.empty:
        for _, g in gene_tss.iterrows():
            tss_lookup[str(g["symbol"]).upper()] = (str(g["chrom"]), int(g["tss"]))

    carollia_aligner = None
    if carollia_seqs is not None and carollia_genome is not None:
        carollia_aligner = _load_carollia_aligner(carollia_genome, carollia_preset)

    carry = [
        c
        for c in ("nearest_gene", "te_name", "te_family", "te_flagged", "nominated")
        if c in candidates.columns
    ]
    rows: list[dict[str, object]] = []
    for _, cand in candidates.iterrows():
        name = cand["name"]
        chrom, start, end = str(cand["chrom"]), int(cand["start"]), int(cand["end"])
        gene = str(cand.get("nearest_gene", "")).upper()
        gene_chrom, gene_pos = tss_lookup.get(gene, (None, None))

        row: dict[str, object] = {"name": name}
        for col in carry:
            row[col] = cand[col]

        for species, lo in lifters.items():
            lift = lift_element(
                lo, species, chrom, start, end, gene_chrom=gene_chrom, gene_pos=gene_pos
            )
            row[f"{species}_class"] = resolve_call(lift)
            row[f"{species}_frac"] = round(lift.fraction, 3)
            row[f"{species}_locus"] = lift.locus
            row[f"{species}_gene_lifts"] = lift.gene_lifts

        if carollia_aligner is not None:
            seq = (carollia_seqs or {}).get(name, "")
            call, ident, cov, locus = _carollia_call(seq, carollia_aligner)
            row["bat_carollia_class"] = call
            row["bat_carollia_ident"] = ident
            row["bat_carollia_cov"] = cov
            row["bat_carollia_locus"] = locus
        rows.append(row)
    return pd.DataFrame(rows)


def convergence_verdict(
    presence: pd.DataFrame,
    *,
    trait_positive: tuple[str, ...] = TRAIT_POSITIVE,
    trait_negative: tuple[str, ...] = TRAIT_NEGATIVE,
) -> pd.DataFrame:
    """Label each element's cross-species pattern against the trait phylogeny.

    Adds a ``verdict``:

    * **CONVERGENT** — PRESENT in every trait-positive lineage *and* lost
      (ABSENT/DEGRADED) in every trait-negative lineage: presence tracks the
      trait, not the tree. (Presence-correlation alone cannot prove
      *independent* gain — shared ancestry then rodent loss yields the same
      pattern — so this is a strong nomination, adjudicated in the report.)
    * **CONSERVED** — PRESENT in trait-positive *and* trait-negative
      lineages: an ancestral element, not a convergent trigger.
    * **ABSENT** — not recoverable in any non-human lineage.
    * **MIXED** — any other pattern (partial / GAP-confounded / ambiguous).

    GAP calls (assembly/synteny gaps) are treated as *not present* but are
    not counted as informative loss, so a GAP in a trait-negative lineage
    cannot by itself manufacture a CONVERGENT verdict.
    """
    df = presence.copy()

    def call(row: pd.Series) -> str:
        pos = [str(row.get(f"{sp}_class", "GAP")) for sp in trait_positive]
        neg = [str(row.get(f"{sp}_class", "GAP")) for sp in trait_negative]
        pos_present = bool(pos) and all(c == "PRESENT" for c in pos)
        neg_lost = bool(neg) and all(c in ("ABSENT", "DEGRADED") for c in neg)
        neg_present = bool(neg) and all(c == "PRESENT" for c in neg)
        all_gone = all(c in ("ABSENT", "GAP") for c in pos + neg)
        if pos_present and neg_lost:
            return "CONVERGENT"
        if pos_present and neg_present:
            return "CONSERVED"
        if all_gone:
            return "ABSENT"
        return "MIXED"

    df["verdict"] = df.apply(call, axis=1)
    order = {"CONVERGENT": 0, "CONSERVED": 1, "MIXED": 2, "ABSENT": 3}
    sort_keys = ["verdict"]
    ascending = [True]
    if "nominated" in df.columns:
        sort_keys.append("nominated")
        ascending.append(False)
    df = df.sort_values(
        sort_keys,
        key=lambda s: s.map(order) if s.name == "verdict" else s,
        ascending=ascending,
    ).reset_index(drop=True)
    return df
