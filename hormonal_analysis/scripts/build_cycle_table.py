"""Normalize seed and raw inputs into one long-form cycle table.

Reads:
- hormonal_analysis/data/seed/cycle_seed.csv
- hormonal_analysis/sources.yaml (for source_id validation)

Writes:
- hormonal_analysis/data/processed/cycle_long.csv
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd
import yaml

HERE = Path(__file__).resolve().parent.parent
SEED = HERE / "data" / "seed" / "cycle_seed.csv"
SOURCES = HERE / "sources.yaml"
OUT = HERE / "data" / "processed" / "cycle_long.csv"

REQUIRED_COLUMNS = [
    "species",
    "source_id",
    "hormone",
    "coordinate_type",
    "coordinate",
    "value",
    "unit",
]
ALL_COLUMNS = REQUIRED_COLUMNS + ["lower", "upper", "n", "notes"]
VALID_COORDINATE_TYPES = {"day", "stage", "normalized"}


def _known_source_ids() -> set[str]:
    with SOURCES.open() as fh:
        manifest = yaml.safe_load(fh)
    ids: set[str] = set()
    for section in ("cycle_sources", "contraception_sources"):
        for entry in manifest.get(section, []) or []:
            ids.add(entry["id"])
    return ids


def main() -> int:
    if not SEED.exists():
        raise SystemExit(f"missing seed file: {SEED}")

    df = pd.read_csv(SEED)
    missing = [c for c in REQUIRED_COLUMNS if c not in df.columns]
    if missing:
        raise SystemExit(f"seed missing required columns: {missing}")

    for col in ("lower", "upper", "n", "notes"):
        if col not in df.columns:
            df[col] = pd.NA

    bad_types = sorted(set(df["coordinate_type"]) - VALID_COORDINATE_TYPES)
    if bad_types:
        raise SystemExit(f"invalid coordinate_type values: {bad_types}")

    known = _known_source_ids()
    unknown_sources = sorted(set(df["source_id"]) - known)
    if unknown_sources:
        raise SystemExit(f"source_id not in sources.yaml: {unknown_sources}")

    out = df[ALL_COLUMNS].copy()
    OUT.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(OUT, index=False)
    print(f"wrote {len(out)} rows -> {OUT.relative_to(HERE.parent)}")
    by_species = out.groupby("species").size().to_dict()
    for species, n in sorted(by_species.items()):
        print(f"  {species}: {n} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
