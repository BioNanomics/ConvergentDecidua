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

# Canonical ordering of mouse and rat estrous stages.
STAGE_ORDER = ["proestrus", "estrus", "metestrus", "diestrus"]

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
        if coord in STAGE_ORDER:
            return STAGE_ORDER.index(coord) / (len(STAGE_ORDER) - 1)
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
        sub = sub[sub["coordinate"].isin(STAGE_ORDER)].copy()
        sub["x"] = sub["coordinate"].map(lambda s: STAGE_ORDER.index(s))
        x_label = "Estrous stage"
        x_ticks = (list(range(len(STAGE_ORDER))), STAGE_ORDER)
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
