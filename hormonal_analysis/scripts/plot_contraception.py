"""Generate the birth-control endogenous-impact plot.

Reads:
- hormonal_analysis/data/seed/contraception_endogenous_seed.csv
- hormonal_analysis/sources.yaml (source_id validation)

Writes:
- hormonal_analysis/data/processed/contraception_endogenous.csv
- hormonal_analysis/plots/contraception_endogenous.png

The plot is a method x hormone heatmap of the ordinal effect_score
(-3 strongly_suppressed .. +3 elevated). This is the Phase 6A
deliverable. Phase 6B (PK curves) is deferred.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

HERE = Path(__file__).resolve().parent.parent
SEED = HERE / "data" / "seed" / "contraception_endogenous_seed.csv"
SOURCES = HERE / "sources.yaml"
PROCESSED = HERE / "data" / "processed" / "contraception_endogenous.csv"
PLOT = HERE / "plots" / "contraception_endogenous.png"

METHOD_ORDER = [
    "combined_oral_contraceptive",
    "vaginal_ring",
    "transdermal_patch",
    "etonogestrel_implant",
    "dmpa",
    "progestin_only_pill",
    "levonorgestrel_iud",
]
HORMONE_ORDER = ["lh", "fsh", "estradiol", "progesterone"]
EFFECT_TICK_LABELS = {
    -3: "strongly suppressed",
    -2: "moderately suppressed",
    -1: "mildly suppressed",
    0: "unchanged",
    1: "mildly elevated",
    2: "moderately elevated",
    3: "strongly elevated",
}


def _known_source_ids() -> set[str]:
    with SOURCES.open() as fh:
        manifest = yaml.safe_load(fh)
    ids: set[str] = set()
    for section in ("cycle_sources", "contraception_sources"):
        for entry in manifest.get(section, []) or []:
            ids.add(entry["id"])
    return ids


def _pretty(method: str) -> str:
    return method.replace("_", " ")


def main() -> int:
    if not SEED.exists():
        raise SystemExit(f"missing seed file: {SEED}")

    df = pd.read_csv(SEED)
    known = _known_source_ids()
    bad = sorted(set(df["source_id"]) - known)
    if bad:
        raise SystemExit(f"source_id not in sources.yaml: {bad}")

    PROCESSED.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(PROCESSED, index=False)
    print(f"wrote {len(df)} rows -> {PROCESSED.relative_to(HERE.parent)}")

    pivot = df.pivot_table(
        index="method", columns="hormone", values="effect_score", aggfunc="mean"
    ).reindex(index=METHOD_ORDER, columns=HORMONE_ORDER)

    fig, ax = plt.subplots(figsize=(7, 5))
    im = ax.imshow(pivot.values, cmap="RdBu", vmin=-3, vmax=3, aspect="auto")
    ax.set_xticks(range(len(HORMONE_ORDER)))
    ax.set_xticklabels([h.upper() for h in HORMONE_ORDER])
    ax.set_yticks(range(len(METHOD_ORDER)))
    ax.set_yticklabels([_pretty(m) for m in METHOD_ORDER])
    ax.set_title("Birth control: impact on endogenous cycle hormones")

    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            value = pivot.values[i, j]
            if np.isnan(value):
                continue
            ax.text(
                j,
                i,
                f"{int(value):+d}",
                ha="center",
                va="center",
                color="black" if abs(value) < 2 else "white",
                fontsize=10,
            )

    cbar = fig.colorbar(im, ax=ax, ticks=sorted(EFFECT_TICK_LABELS))
    cbar.ax.set_yticklabels([EFFECT_TICK_LABELS[v] for v in sorted(EFFECT_TICK_LABELS)])
    cbar.set_label("Effect on endogenous hormone")

    fig.tight_layout()
    PLOT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOT, dpi=150)
    plt.close(fig)
    print(f"wrote {PLOT.relative_to(HERE.parent)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
