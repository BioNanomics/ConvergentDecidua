# Long-form cycle hormone schema

Single processed table at `data/processed/cycle_long.csv`. One row per
(species, study, hormone, cycle coordinate). All seed and raw inputs
are normalized into this shape.

## Columns

| column | type | meaning |
|---|---|---|
| `species` | str | Matches a row in `species_matrix.csv`. |
| `source_id` | str | Matches an `id` in `sources.yaml` `cycle_sources`. |
| `hormone` | str | One of `estradiol`, `progesterone`, `lh`, `fsh`, `prolactin`. |
| `coordinate_type` | str | `day` (1-indexed cycle day), `stage` (named phase), or `normalized` (0–1 cycle fraction). |
| `coordinate` | str | Numeric for `day` / `normalized`, name for `stage` (e.g. `proestrus`, `mid_luteal`). |
| `value` | float | Reported central tendency (mean or median). |
| `unit` | str | Original assay unit (e.g. `pg/mL`, `ng/mL`, `mIU/mL`, `nmol/L`). Preserved, **not** silently rescaled. |
| `lower` | float | Optional lower bound (SE, SD, or range low). May be empty. |
| `upper` | float | Optional upper bound. May be empty. |
| `n` | int | Subject count if reported. May be empty. |
| `notes` | str | Free-form, e.g. assay method or population caveats. |

## Rules

- **Never convert units silently.** If a comparison needs matched
  units, do the conversion in the plot script and label it.
- **Preserve coordinate type.** A mouse stage axis and a human day
  axis are not interchangeable; cross-species comparison plots use
  the `normalized` coordinate explicitly.
- **Every row must have `source_id`** that exists in `sources.yaml`.
  `build_cycle_table.py` enforces this.
