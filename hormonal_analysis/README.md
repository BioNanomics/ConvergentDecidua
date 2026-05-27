# hormonal_analysis

Comparative reproductive-cycle hormone data, cross-species cycle
plots, and contraceptive-method cycle-impact plots.

See **[results.md](results.md)** for the rendered figures and
discussion.

## What's here

- **[species_matrix.csv](species_matrix.csv)** — compact matrix
  (phylogenetic closeness, spontaneous decidualization, menstruation,
  genome quality, hormone-data quality, downloadability, status).
  Rendered version with per-species framing notes is
  [results.md §1](results.md#1-species-matrix-at-a-glance).
- **[sources.yaml](sources.yaml)** — every value's provenance, with an
  `access_class` of `automated_download`, `manual_table_extraction`,
  or `reference_only`.
- **[schema/cycle_long.md](schema/cycle_long.md)** /
  **[schema/contraception.md](schema/contraception.md)** — column
  definitions for the processed tables.
- **data/seed/** — hand-curated, citation-tagged seed CSVs (the
  current source of truth; see caveats in `results.md`).
- **data/processed/** — normalized long-form tables built from seed.
- **plots/** — generated PNGs, committed via Git LFS.
- **scripts/** — standalone Python (pandas, pyyaml, matplotlib only).

## Layout

```
hormonal_analysis/
  README.md
  results.md
  species_matrix.csv
  sources.yaml
  schema/
    cycle_long.md
    contraception.md
  data/
    raw/         downloaded source files (gitignored)
    processed/   normalized long-form tables (gitignored, regenerable)
    seed/        hand-curated seed CSVs (tracked)
  plots/         generated PNGs (tracked via Git LFS)
  scripts/
    fetch_sources.py
    build_cycle_table.py
    plot_cycles.py
    plot_contraception.py
```

## Reproduce

From the repo root, with the project venv active:

```bash
python hormonal_analysis/scripts/fetch_sources.py        # no-op until automated sources are added
python hormonal_analysis/scripts/build_cycle_table.py    # seed -> data/processed/cycle_long.csv
python hormonal_analysis/scripts/plot_cycles.py          # per-species + cross-species cycle plots
python hormonal_analysis/scripts/plot_contraception.py   # contraception endogenous-impact heatmap
```

The scripts are deterministic and idempotent. Outputs land in
`data/processed/` and `plots/`.
