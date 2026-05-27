"""Generate cycle hormone plots.

Produces:
- hormonal_analysis/plots/cycle_<species>.png    (one per plotted species)
- hormonal_analysis/plots/cycle_cross_species_<hormone>.png
  (cross-species comparison on a normalized 0..1 cycle axis)

Cross-species comparisons project each species onto a normalized
cycle position so menstrual and estrous cycles can be aligned without
pretending they share day semantics. Units are preserved per row;
plot legends label the unit.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

HERE = Path(__file__).resolve().parent.parent
LONG = HERE / "data" / "processed" / "cycle_long.csv"
PLOTS = HERE / "plots"

# Canonical ordering of mouse and rat estrous stages (used for the
# per-species plot, where the x-axis is by convention proestrus -> estrus
# -> metestrus -> diestrus).
STAGE_ORDER = ["proestrus", "estrus", "metestrus", "diestrus"]

# Display order for the per-species rodent plot. We DEPART from the
# textbook proestrus -> diestrus order so the rodent plots can be
# visually compared panel-to-panel against the human cycle plot in
# section 2.1 of results.md. With this ordering, proestrus (the
# rodent pre-ovulatory LH surge) sits in the middle of the plot at
# the same visual position as the human day-13 LH surge, instead of
# at the left edge. Readers should be aware this is NOT the
# textbook estrous-cycle layout -- it is a deliberate alignment
# choice for cross-panel readability and is annotated on the plot
# itself.
STAGE_DISPLAY_ORDER = ["diestrus", "proestrus", "estrus", "metestrus"]

# Canonical rodent estrous-cycle stage durations, in hours. Values
# are the Marcondes/Bianchi/Tanno 2002 rat 4-day cycle (Braz J Biol
# 62(4A):609-614, Table I); mouse durations are similar within
# ~10-20 % so we treat these as a common rodent default rather than
# curating species-specific durations (which would require carving
# the seed into separate mouse_estrous_canonical_durations and
# rat_estrous_canonical_durations). Total = 102 h ~= 4.25 d.
STAGE_DURATIONS_H = {
    "proestrus": 12,
    "estrus": 12,
    "metestrus": 21,
    "diestrus": 57,
}

# Biological cycle order (diestrus -> proestrus -> estrus -> metestrus
# -> diestrus). We build duration-weighted stage midpoints on the
# 0..1 cycle, then rotate the whole cycle so that the proestrus
# midpoint sits at 0.43 -- matching the human day-13 LH surge at
# (13 - 1) / (28 - 1) = 0.444. The resulting positions are:
#   diestrus  ~ 0.09  (long: 0.559 of cycle wall-clock)
#   proestrus ~ 0.43  (short: 0.118)
#   estrus    ~ 0.55  (short: 0.118)
#   metestrus ~ 0.71  (medium: 0.206)
# Proestrus and estrus are visually tight in cross-species panels
# because the rodent surge + ovulation window is only ~24 h of a
# ~4-day cycle; diestrus correspondingly stretches over more than
# half the normalized axis. This is intentional -- the goal of the
# duration weighting is to make cross-species panels reflect
# biological wall-clock time, not stage count. The PER-SPECIES
# rodent panels in section 2.2 / 2.3 still use equal-width stages
# (one tick per stage); the duration-weighted positions are applied
# only to the cross-species panels in section 3.
_CYCLE_ORDER = ["diestrus", "proestrus", "estrus", "metestrus"]
_TOTAL_H = sum(STAGE_DURATIONS_H.values())
_PROESTRUS_ANCHOR = 0.43


def _compute_stage_normalized() -> dict[str, float]:
    midpoints: dict[str, float] = {}
    cursor = 0.0
    for stage in _CYCLE_ORDER:
        width = STAGE_DURATIONS_H[stage] / _TOTAL_H
        midpoints[stage] = cursor + width / 2
        cursor += width
    shift = _PROESTRUS_ANCHOR - midpoints["proestrus"]
    return {s: round((m + shift) % 1.0, 4) for s, m in midpoints.items()}


STAGE_NORMALIZED = _compute_stage_normalized()

# Human cycle length used to normalize day -> 0..1.
HUMAN_CYCLE_LENGTH_DAYS = 28

HORMONE_LABELS = {
    "estradiol": "Estradiol",
    "progesterone": "Progesterone",
    "lh": "LH",
    "fsh": "FSH",
    "prolactin": "Prolactin",
}


def _normalize_position(row: pd.Series) -> float | None:
    """Map a row's coordinate to a 0..1 cycle position."""
    coord_type = row["coordinate_type"]
    coord = row["coordinate"]
    if coord_type == "day":
        try:
            day = float(coord)
        except (TypeError, ValueError):
            return None
        return max(0.0, min(1.0, (day - 1) / (HUMAN_CYCLE_LENGTH_DAYS - 1)))
    if coord_type == "stage":
        if coord in STAGE_NORMALIZED:
            return STAGE_NORMALIZED[coord]
        return None
    if coord_type == "normalized":
        try:
            return float(coord)
        except (TypeError, ValueError):
            return None
    return None


def _plot_species(species: str, df: pd.DataFrame) -> Path:
    sub = df[df["species"] == species].copy()
    coord_type = sub["coordinate_type"].iloc[0]

    if coord_type == "day":
        sub["x"] = sub["coordinate"].astype(float)
        x_label = "Cycle day"
        x_ticks = None
    elif coord_type == "stage":
        sub = sub[sub["coordinate"].isin(STAGE_DISPLAY_ORDER)].copy()
        sub["x"] = sub["coordinate"].map(lambda s: STAGE_DISPLAY_ORDER.index(s))
        x_label = (
            "Estrous stage (non-textbook order: diestrus first, so the\n"
            "proestrus LH surge aligns with the human day-13 surge in §2.1)"
        )
        x_ticks = (list(range(len(STAGE_DISPLAY_ORDER))), STAGE_DISPLAY_ORDER)
    elif coord_type == "normalized":
        sub["x"] = sub["coordinate"].astype(float)
        x_label = "Normalized cycle position (0 = start, 1 = end)"
        x_ticks = None
    else:
        raise ValueError(f"unsupported coordinate_type for {species}: {coord_type}")

    fig, ax = plt.subplots(figsize=(8, 5))
    for hormone, hdf in sub.groupby("hormone"):
        hdf = hdf.sort_values("x")
        unit = hdf["unit"].iloc[0]
        label = f"{HORMONE_LABELS.get(hormone, hormone)} ({unit})"
        ax.plot(hdf["x"], hdf["value"], marker="o", label=label)
    ax.set_title(f"Reproductive cycle hormones — {species}")
    ax.set_xlabel(x_label)
    ax.set_ylabel("Concentration (see legend for units)")
    if x_ticks:
        ax.set_xticks(x_ticks[0])
        ax.set_xticklabels(x_ticks[1])
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()

    out = PLOTS / f"cycle_{species}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def _plot_cross_species(hormone: str, df: pd.DataFrame) -> Path | None:
    sub = df[df["hormone"] == hormone].copy()
    sub["pos"] = sub.apply(_normalize_position, axis=1)
    sub = sub.dropna(subset=["pos"])
    if sub.empty:
        return None

    fig, ax = plt.subplots(figsize=(8, 5))
    for species, sdf in sub.groupby("species"):
        sdf = sdf.sort_values("pos")
        unit = sdf["unit"].iloc[0]
        label = f"{species} ({unit})"
        ax.plot(sdf["pos"], sdf["value"], marker="o", label=label)

    ax.set_title(f"{HORMONE_LABELS.get(hormone, hormone)} — cross-species (normalized cycle)")
    ax.set_xlabel("Normalized cycle position (0 = start, 1 = end)")
    ax.set_ylabel("Concentration (units vary; see legend)")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best", fontsize=9)
    fig.tight_layout()
    out = PLOTS / f"cycle_cross_species_{hormone}.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    return out


def main() -> int:
    if not LONG.exists():
        raise SystemExit(f"missing processed table: {LONG} (run build_cycle_table.py first)")

    df = pd.read_csv(LONG)
    PLOTS.mkdir(parents=True, exist_ok=True)

    species_list = sorted(df["species"].unique())
    print(f"per-species plots ({len(species_list)}):")
    for species in species_list:
        out = _plot_species(species, df)
        print(f"  {out.relative_to(HERE.parent)}")

    hormones = sorted(df["hormone"].unique())
    print(f"cross-species plots ({len(hormones)}):")
    for hormone in hormones:
        out = _plot_cross_species(hormone, df)
        if out is None:
            print(f"  (skipped) {hormone}: no normalizable rows")
        else:
            print(f"  {out.relative_to(HERE.parent)}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
