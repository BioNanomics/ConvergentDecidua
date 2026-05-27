# Contraception schema

Two tables under `data/processed/`.

## `contraception_endogenous.csv` (Phase 6A)

How each method shifts endogenous LH / FSH / estradiol / progesterone
versus a natural cycle baseline. Qualitative-but-defensible: values
are typically reported in labels and review summaries as percentages
of baseline or as "suppressed below detection."

| column | type | meaning |
|---|---|---|
| `method` | str | e.g. `combined_oral_contraceptive`, `levonorgestrel_iud`, `etonogestrel_implant`, `dmpa`, `vaginal_ring`, `transdermal_patch`, `progestin_only_pill`. |
| `hormone` | str | `lh`, `fsh`, `estradiol`, `progesterone`. |
| `effect` | str | One of `strongly_suppressed`, `moderately_suppressed`, `mildly_suppressed`, `unchanged`, `elevated`. |
| `effect_score` | int | -3..+3 ordinal encoding of `effect` for plotting. |
| `source_id` | str | Matches `sources.yaml` `contraception_sources`. |
| `notes` | str | Free-form. |

## `contraception_pk.csv` (Phase 6B, deferred)

Pharmacokinetic summary per active compound and formulation.

| column | type | meaning |
|---|---|---|
| `method` | str | Same vocabulary as above. |
| `compound` | str | e.g. `ethinylestradiol`, `levonorgestrel`, `etonogestrel`, `medroxyprogesterone_acetate`. |
| `cmax` | float | Peak concentration. |
| `cmax_unit` | str | e.g. `pg/mL`. |
| `tmax_hours` | float | Time to Cmax. |
| `half_life_hours` | float | Elimination half-life. |
| `auc` | float | Optional AUC. |
| `auc_unit` | str | Optional. |
| `source_id` | str | Matches `sources.yaml`. |
| `notes` | str | Free-form. |

Phase 6B is intentionally not seeded yet. See `README.md` "Scope
boundary."
