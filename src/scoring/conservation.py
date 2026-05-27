"""Conserved vs divergent module classification (Q3.3).

Consumes the permutation-FDR table produced by ``src/scoring/null.py``
and assigns each (module × group) a class based on the cross-species
significance pattern.

Class labels (per (module, group) row):

- ``conserved-up`` — both species FDR < threshold and observed_mean -
  null_mean > 0 in both.
- ``conserved-down`` — both species FDR < threshold and observed_mean -
  null_mean < 0 in both.
- ``divergent`` — both species FDR < threshold but with opposite signs.
- ``human-biased-up`` / ``human-biased-down`` — human FDR < threshold
  only (with the matching sign); mouse not significant.
- ``mouse-biased-up`` / ``mouse-biased-down`` — symmetric.
- ``neutral`` — neither species significant.

The module-level summary in :func:`summarise_modules` rolls up to one
row per module by reporting the strongest pattern across groups (max
abs deviation among rows tagged conserved or biased).
"""

from __future__ import annotations

import pandas as pd

CONSERVATION_COLS = ("module", "species", "group", "observed_mean", "null_mean", "fdr")


def classify_conservation(fdr_table: pd.DataFrame, fdr_threshold: float = 0.05) -> pd.DataFrame:
    """Pivot the long FDR table to one row per (module, group) with a class."""
    missing = [c for c in CONSERVATION_COLS if c not in fdr_table.columns]
    if missing:
        msg = f"fdr_table missing columns: {missing}"
        raise ValueError(msg)

    work = fdr_table.copy()
    work["deviation"] = work["observed_mean"] - work["null_mean"]
    work["sig"] = work["fdr"] < fdr_threshold
    work["direction"] = work["deviation"].apply(lambda d: "up" if d > 0 else "down")

    species_levels = sorted(work["species"].unique())
    rows: list[dict[str, object]] = []
    for (module, group), grp in work.groupby(["module", "group"]):
        per_species = grp.set_index("species")
        per_species = per_species.reindex(species_levels)

        sig_flags = {
            sp: bool(per_species.loc[sp, "sig"]) if sp in grp["species"].values else False
            for sp in species_levels
        }
        dirs = {
            sp: per_species.loc[sp, "direction"] if sp in grp["species"].values else None
            for sp in species_levels
        }
        devs = {
            sp: float(per_species.loc[sp, "deviation"]) if sp in grp["species"].values else 0.0
            for sp in species_levels
        }

        sig_species = [sp for sp, s in sig_flags.items() if s]
        if len(sig_species) == len(species_levels) and len(species_levels) >= 2:
            if len({dirs[sp] for sp in sig_species}) == 1:
                cls = f"conserved-{dirs[sig_species[0]]}"
            else:
                cls = "divergent"
        elif len(sig_species) == 1:
            sp = sig_species[0]
            cls = f"{sp}-biased-{dirs[sp]}"
        else:
            cls = "neutral"

        row: dict[str, object] = {"module": module, "group": group, "class": cls}
        for sp in species_levels:
            row[f"{sp}_deviation"] = devs[sp]
            row[f"{sp}_fdr"] = (
                float(per_species.loc[sp, "fdr"]) if sp in grp["species"].values else float("nan")
            )
        rows.append(row)

    return pd.DataFrame(rows).sort_values(["module", "group"]).reset_index(drop=True)


def summarise_modules(conservation_table: pd.DataFrame) -> pd.DataFrame:
    """One row per module: strongest pattern across all groups.

    A module is labelled ``conserved-up`` / ``conserved-down`` if any
    group is conserved in that direction. Otherwise it falls back to
    the species-biased class with the largest |deviation| in either
    species, or ``neutral`` if no group reached significance.
    """
    rows: list[dict[str, object]] = []
    for module, grp in conservation_table.groupby("module"):
        conserved_rows = grp[grp["class"].str.startswith("conserved-")]
        if not conserved_rows.empty:
            # Pick the conserved direction with the strongest mean deviation.
            best = conserved_rows.copy()
            dev_cols = [c for c in best.columns if c.endswith("_deviation")]
            best["max_abs_dev"] = best[dev_cols].abs().max(axis=1)
            top = best.sort_values("max_abs_dev", ascending=False).iloc[0]
            rows.append(
                {
                    "module": module,
                    "summary_class": top["class"],
                    "top_group": top["group"],
                    "max_abs_deviation": float(top["max_abs_dev"]),
                    "n_conserved_groups": int(len(conserved_rows)),
                }
            )
            continue

        biased = grp[grp["class"].str.contains("-biased-")]
        if not biased.empty:
            dev_cols = [c for c in biased.columns if c.endswith("_deviation")]
            biased = biased.copy()
            biased["max_abs_dev"] = biased[dev_cols].abs().max(axis=1)
            top = biased.sort_values("max_abs_dev", ascending=False).iloc[0]
            rows.append(
                {
                    "module": module,
                    "summary_class": top["class"],
                    "top_group": top["group"],
                    "max_abs_deviation": float(top["max_abs_dev"]),
                    "n_conserved_groups": 0,
                }
            )
            continue

        rows.append(
            {
                "module": module,
                "summary_class": "neutral",
                "top_group": "",
                "max_abs_deviation": 0.0,
                "n_conserved_groups": 0,
            }
        )

    return pd.DataFrame(rows).sort_values("module").reset_index(drop=True)
