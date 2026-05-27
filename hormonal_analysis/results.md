# Results — comparative cycle hormones + contraception impact

> **Status:** first-pass results from seeded canonical literature
> values, not from raw downloaded datasets. All quantitative claims
> are at the **shape-of-curve** level. See
> [README.md](README.md#scope-boundary) for the scope boundary that
> separates this workspace from the year-one atlas core.

## 1. Species matrix at a glance

The reproducible core is **human + mouse + rat** (trait-positive
human; trait-negative mouse, rat). Other rows in
[species_matrix.csv](species_matrix.csv) are exploratory comparators
and carry caveats in [species_matrix.md](species_matrix.md):

- **Rhesus macaque** — menstruating, spontaneously decidualizing. On
  the *wrong* axis for a trait-negative control; included as
  `reference_only` for future menstruation-but-induced contrasts.
- **Common marmoset** — `unconfirmed_likely_no` for spontaneous
  decidualization; the literature is ambiguous and we have not yet
  pinned down a definitive citation either way.
- **Gray mouse lemur** — cleaner non-menstruating, non-spontaneous
  reference than mouse phylogenetically, but hormone time-series data
  is sparse; flagged as exploratory.
- **Spiny mouse** — *menstruates* but is muroid; the most
  evolutionarily interesting comparator and the row most likely to
  graduate from exploratory to core if a usable hormone series can be
  curated.
- **Bat** — intentionally absent; out-of-scope per
  [PLAN.md](../PLAN.md) locked decisions.

## 2. Cycle hormone curves — per species, native axis

Each per-species plot is on its **native time coordinate**: human on
days 1–28, mouse and rat on the four-stage estrous axis
(proestrus → estrus → metestrus → diestrus). Units are preserved per
hormone (E2 in pg/mL; P4 in ng/mL; gonadotropins in mIU/mL for human
and ng/mL for rodents — **not** unit-converted across species).

### 2.1 Human (28-day cycle)

![Human menstrual cycle hormones](plots/cycle_human.png)

- E2 shows the canonical two-peak shape: a sharp pre-ovulatory peak
  around day 12 (~300 pg/mL) and a broader mid-luteal secondary peak
  around day 21 (~180 pg/mL).
- P4 stays at baseline (<1 ng/mL) through the entire follicular phase
  and only rises after ovulation, peaking mid-luteal (~12 ng/mL) and
  collapsing premenstrually.
- The LH surge (~55 mIU/mL on day 13) is sharp and brief; the FSH
  peak co-occurs but is smaller (~18 mIU/mL).

### 2.2 Mouse (4–5 day estrous cycle)

![Mouse estrous cycle hormones](plots/cycle_mouse.png)

- E2 peaks at proestrus (~40 pg/mL) and falls through estrus into
  metestrus.
- P4 peaks at **metestrus** (~15 ng/mL) — not at proestrus — reflecting
  the brief, induced (not spontaneous) corpus luteum.
- LH and FSH surges are both concentrated at proestrus.

### 2.3 Rat (4–5 day estrous cycle)

![Rat estrous cycle hormones](plots/cycle_rat.png)

- Shape is broadly similar to mouse but with a higher proestrus P4
  spike (~30 ng/mL) reflecting the rat's stronger pre-ovulatory P4.
- The defining rat feature included here is the **prolactin surge at
  proestrus** (~80 ng/mL) — paired with the LH surge and absent from
  the mouse seed. This is the one rat-specific axis curated in the
  seed.

## 3. Cross-species comparison — normalized cycle position

Cross-species plots place each species on a common 0..1 cycle-position
axis (human: day / 28; rodents: stage index / 3). **Units are still
per-species** so absolute magnitudes are not directly comparable —
these plots are for **phase alignment**, not amplitude comparison.

### 3.1 Estradiol

![Estradiol across species](plots/cycle_cross_species_estradiol.png)

The pre-ovulatory E2 peak lands at roughly the same normalized
position across all three species (~40–45 % through cycle), supporting
phase alignment but **not** evidence of magnitude equivalence (human
peak is 6–7× higher than rodent in absolute pg/mL).

### 3.2 Progesterone

![Progesterone across species](plots/cycle_cross_species_progesterone.png)

This is the most biologically informative cross-species panel: the
human luteal P4 plateau is **broad** (days 14–25, roughly 50–90 % of
cycle), whereas rodent P4 is a brief metestrus / post-ovulatory
bump. The width of the luteal P4 window — not its peak height — is
the candidate correlate of spontaneous decidualization. **This is a
hypothesis-generating observation, not a tested claim.**

### 3.3 LH and FSH

![LH across species](plots/cycle_cross_species_lh.png)
![FSH across species](plots/cycle_cross_species_fsh.png)

LH and FSH surges co-localize at the same normalized cycle position
(~40–45 %) in all three species. Shape is consistent; this is the
expected cross-species conservation and serves mainly as a sanity
check that the normalization is sensible.

### 3.4 Prolactin

![Prolactin across species](plots/cycle_cross_species_prolactin.png)

Only the rat seed currently includes a prolactin time series.
Recovering a comparable human and mouse trace is a clear next step
before any prolactin-based cross-species claim.

## 4. Contraception → endogenous-hormone impact

![Contraception endogenous hormone impact](plots/contraception_endogenous.png)

Ordinal **effect_score ∈ {−3, −2, −1, 0, +1, +2, +3}** per method ×
endogenous hormone, sourced from FDA drug labels
([sources.yaml](sources.yaml) `fda_drug_labels`). Negative = the
method *suppresses the endogenous hormone*; this does **not** track
exogenous hormone exposure (e.g., ethinyl estradiol from a COC is
*not* counted as "endogenous E2").

Patterns the figure makes obvious:

- **Combined oral, vaginal ring, transdermal patch** form a cluster —
  endogenous LH/FSH/P4 strongly suppressed (−3), endogenous E2
  moderately suppressed (−2). They share an HPO-axis-suppression
  mechanism.
- **DMPA and the etonogestrel implant** also suppress endogenous P4
  strongly via reliable ovulation suppression, but the implant
  preserves endogenous E2 better than DMPA.
- **Levonorgestrel IUD** sits at 0 across the row — ovulatory cycles
  are largely preserved; its mechanism is endometrial / cervical, not
  HPO-axis.
- **Progestin-only pill** is the intermediate row — variable
  ovulation suppression yields mostly mild (−1) endogenous effects.

The figure is therefore primarily a **mechanism-of-action
classifier**, not a clinical efficacy comparison.

## 5. Caveats — read before citing any number

- **Seed values are canonical literature estimates** transcribed by
  hand. They reflect textbook / review-article means and are
  appropriate for shape-of-curve discussion only. They are **not**
  individual-level or even study-level data. Real data ingestion
  requires manual table extraction from the references listed in
  [sources.yaml](sources.yaml).
- **Units are not silently converted across species** (gonadotropin
  units in particular differ; do not read absolute amplitudes off
  cross-species plots).
- **Cycle position normalization is uniform-spacing on stage index for
  rodents**, which compresses the brief estrus phase and stretches
  the longer diestrus phase relative to wall-clock time. A
  duration-weighted normalization is a sensible next iteration.
- **The contraception matrix is endogenous-hormone-only.** Exogenous
  PK curves (ethinyl estradiol levels, levonorgestrel serum
  concentration vs time) are Phase 6B and deliberately deferred.
- **Primate and spiny mouse rows are exploratory.** Do not let them
  carry quantitative claims.

## 6. Suggested next iterations (none committed)

In rough order of yield-per-effort:

1. **Spiny mouse cycle series** — if a usable hormone time series can
   be curated from Bellofiore et al. and follow-up work, this becomes
   the most informative single addition: a menstruating muroid is the
   one comparator that can separate "muroid lineage" from "spontaneous
   decidualization" as the explanatory axis.
2. **Duration-weighted rodent cycle normalization** — replace the
   stage-index normalization with stage-duration weights (estrus is
   ~12 h, diestrus ~48 h) so cross-species panels reflect biological
   time, not stage count.
3. **Human prolactin series** to complete the cross-species prolactin
   panel.
4. **Contraception Phase 6B (exogenous PK curves)** — only if
   regulatory or evolutionary questions explicitly need it; the
   endogenous-impact matrix is sufficient for atlas-context framing
   alone.

None of these are commitments. Each would be a new scoped task,
proposed against PLAN.md, not a silent expansion of the year-one
atlas.
