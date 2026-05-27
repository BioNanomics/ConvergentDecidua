# Results — comparative cycle hormones + contraception impact

> **Status:** first-pass results from seeded canonical literature
> values, not from raw downloaded datasets. All quantitative claims
> are at the **shape-of-curve** level. See
> [README.md](README.md#scope-boundary) for the scope boundary that
> separates this workspace from the year-one atlas core.

## 1. Species matrix at a glance

The reproducible core is **human + mouse + rat + spiny mouse**
(trait-positive: human, spiny mouse; trait-negative: mouse, rat).
Other rows in [species_matrix.csv](species_matrix.csv) are
exploratory comparators (see §1.2 for caveats).

### 1.1 Comparison matrix

Compact comparison for the "closely related but lacking spontaneous
decidualization" framing. Columns:

- **phylogenetic_closeness_to_human** — qualitative grouping; not a
  divergence time.
- **spontaneous_decidualization** — `yes` / `no` /
  `unconfirmed_likely_no` / `likely_no`. The strict
  trait-positive / trait-negative axis for this project.
- **menstruation** — related trait but **not** equivalent to
  spontaneous decidualization (see PMID 36304046).
- **genome_quality** — coarse bin from NCBI Datasets at time of
  writing.
- **hormone_data_quality** — practical availability of
  cycle-resolved serum/plasma concentration data.
- **downloadability** — whether usable data is currently behind a
  fetchable URL with a stable schema (`structured_table`), available
  only as a published table or label needing one-time transcription
  (`manual_extraction`), or not usable in the first pass
  (`reference_only`).
- **first_pass_status** — `plotted` species drive the cycle figures;
  `exploratory` species are present in the matrix and source
  manifest but only plotted if a downloadable source is confirmed;
  `reference_only` species inform the comparative framing but are
  not plotted.

| species | closeness | spont. decidu. | menstruation | genome | hormone data | downloadability | first-pass |
|---|---|---|---|---|---|---|---|
| human | reference | yes | yes | T2T complete | excellent, daily serum | structured table | plotted |
| rhesus macaque | close Old World monkey | yes | yes | T2T complete | good, cycle resolved | manual extraction | reference only |
| common marmoset | close New World monkey | unconfirmed, likely no | no | chromosome scale | moderate, ovarian cycle | manual extraction | exploratory |
| gray mouse lemur | strepsirrhine | likely no | no | chromosome scale | sparse | reference only | exploratory |
| spiny mouse (*Acomys*) | rodent (Acomys) | yes | yes | chromosome scale | sparse, phase resolved | manual extraction | plotted |
| mouse | rodent outgroup | no | no | T2T complete | good, stage resolved | manual extraction | plotted |
| rat | rodent outgroup | no | no | chromosome scale | good, stage resolved | manual extraction | plotted |

### 1.2 Per-species framing notes

- **Rhesus macaque** is intentionally **reference only** because it
  sits on the wrong side of the spontaneous-decidualization axis for
  the trait-negative control role, even though it has strong cycle
  hormone data. It is the right comparator for menstruation work but
  the wrong comparator for "closely related, without spontaneous
  decidualization."
- **Common marmoset** is marked **unconfirmed, likely no**: its
  endometrial literature describes cycle-dependent uterine priming
  around implantation without a clean cyclical spontaneous decidual
  program. This is the residual question from the prior research
  turn.
- **Gray mouse lemur** is the cleaner candidate for a phenotype-clean
  primate trait-negative control on the spontaneous-decidualization
  axis, but endocrine data is too sparse for first-pass plotting.
- **Spiny mouse** *menstruates* and spontaneously decidualizes inside
  the muroid clade. Plotted for **P4, prolactin, and E2** using
  Bellofiore et al. 2017 (AJOG) / 2018 (JME) Tables 2-3 for P4 and
  prolactin, plus Bellofiore et al. 2021 (Human Reproduction) Fig 1e
  for estradiol (digitized phase means: ~72 pg/mL menses, ~116 pg/mL
  proliferative peak, ~107 pg/mL secretory secondary elevation; n=9
  6-month females). **LH and FSH absolute serum values cannot be
  measured in this species at present** — no validated cross-
  reactive immunoassay exists and skin autotomy prevents the serial
  blood draws that would be needed for cycle-resolved gonadotropin
  tracking (McKenna 2020 PLoS ONE; McKenna 2021 Sci Rep). This is
  recorded as a species-level technical ceiling in
  [sources.yaml](sources.yaml) (`spiny_mouse_gonadotropin_ceiling`),
  not a missing-paper gap.
- **Bat** is intentionally absent from this matrix because the
  project scope downgraded it to sequence-first per
  [PLAN.md](../PLAN.md) locked decisions; revisit if a downloadable
  bat hormone-cycle table appears.

See [sources.yaml](sources.yaml) for source-level provenance behind
each `hormone data` rating.

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

> **Note on x-axis ordering (non-textbook).** Textbook estrous-cycle
> diagrams conventionally start at **proestrus** and end at
> **diestrus**. The plot above departs from that convention and uses
> the order `diestrus → proestrus → estrus → metestrus`. This is
> a deliberate choice, *not* the biological default, made so that
> the rodent pre-ovulatory LH surge (proestrus) sits in the middle
> of the panel at the same visual position as the human day-13 LH
> surge in §2.1, making the four §2 panels readable as four
> versions of the same phase shape. The cross-species axis in §3
> already does this alignment quantitatively; §2 panels are
> presented this way only for visual consistency. Readers used to
> the textbook convention should mentally rotate the panel one step
> to the left. Stage *durations* are still equal-width (the
> duration-weighted normalization is the §6 next-iteration task).

- E2 peaks at proestrus (~40 pg/mL) and falls through estrus into
  metestrus.
- P4 peaks at **metestrus** (~15 ng/mL) — not at proestrus — reflecting
  the brief, induced (not spontaneous) corpus luteum.
- LH and FSH surges are both concentrated at proestrus.

### 2.3 Rat (4–5 day estrous cycle)

![Rat estrous cycle hormones](plots/cycle_rat.png)

> **Same non-textbook x-axis ordering as §2.2.** The stages are
> drawn as `diestrus → proestrus → estrus → metestrus` to keep
> proestrus centred. See the note under §2.2 for the rationale.

- Shape is broadly similar to mouse but with a higher proestrus P4
  spike (~30 ng/mL) reflecting the rat's stronger pre-ovulatory P4.
- The defining rat feature included here is the **prolactin surge at
  proestrus** (~80 ng/mL) — paired with the LH surge and absent from
  the mouse seed. This is the one rat-specific axis curated in the
  seed.

### 2.4 Spiny mouse (*Acomys cahirinus*, 6–10 day menstrual cycle)

![Spiny mouse cycle hormones](plots/cycle_spiny_mouse.png)

The spiny mouse runs a 6–10 day **menstrual** (not estrous) cycle
with ~3 days of overt menses; phases here are positioned on a
normalized cycle axis because the proestrus/estrus/metestrus/diestrus
labels do not apply.

- **P4** (ng/mL, ELISA, cardiac puncture; Bellofiore 2018 Table 2
  tabulating Bellofiore 2017): early follicular midpoint **~47**
  (range 30–64) → late luteal midpoint **~135** (range 70–199). The
  ~3.1× luteal/follicular fold-increase is exactly at the threshold
  Bellofiore et al. 2018 propose as the unique endocrine feature of
  menstruating species (≥3-fold rise; humans ~5×, rhesus / baboon
  ~6×, fulvous fruit bat ~3×; non-menstruating mouse/rat/sheep <3×
  or no rise).
- **Prolactin** (IU/L, DXI immunoassay): a secondary surge in the
  late luteal phase (~11.7 IU/L) above the early follicular baseline
  (~3.3 IU/L), the only phase difference reported as significant in
  the source (p < 0.05). The spiny mouse pattern mirrors the late-
  luteal prolactin rise seen in humans and rhesus; the mouse and
  sheep prolactin profiles in the same table do not show this rise.
- **Estradiol** (pg/mL, Calbiotech Mouse/Rat ELISA ES180S-100,
  cardiac puncture; Bellofiore 2021 Hum Reprod Fig 1e and Suppl
  Table SI assay validation, n = 9 6-month females). Three phase
  means digitized by eye from Fig 1e (which colours individual
  data points by cycle phase): **menses ~72 pg/mL** (red dots, n=2),
  **proliferative ~116 pg/mL** (orange, n=3), **secretory ~107
  pg/mL** (lavender, n=4). The paper's quoted ~50–140 pg/mL range
  is the within-cohort span across individuals. This is
  qualitatively a human-like two-peak shape (proliferative peak
  with secretory secondary elevation), not the single proestrus
  peak of mouse or rat.
- **Unit note:** spiny mouse prolactin is reported in **IU/L** by
  DXI immunoassay (Bellofiore 2018 Table 3 footnote *), distinct
  from the ng/mL RIA convention used for the rat trace in §3.4. Do
  not read absolute prolactin amplitudes across these two species.
- **LH / FSH — not measurable in this species at present.** Two
  Bellofiore / McKenna group papers state this explicitly: no
  validated spiny mouse LH or FSH immunoassay exists, cross-
  reactive antibodies from other rodents and humans have not
  worked, and skin autotomy precludes the serial blood draws that
  would be needed for cycle-resolved gonadotropin tracking
  (McKenna 2020 PLoS ONE; McKenna 2021 Sci Rep). The gap is a
  species-level technical ceiling, not a missing-paper gap, and is
  recorded in [sources.yaml](sources.yaml) as
  `spiny_mouse_gonadotropin_ceiling`. Closing it requires new
  laboratory assay development, not further reading.

## 3. Cross-species comparison — normalized cycle position

Cross-species plots place each species on a common 0..1 cycle-position
axis. Mapping rules:

- **Human:** `(day - 1) / 27`, so menses-onset = 0.00 and day 28 =
  1.00.
- **Mouse / rat:** stages are placed at **duration-weighted midpoints**
  of their canonical wall-clock durations (proestrus 12 h, estrus 12 h,
  metestrus 21 h, diestrus 57 h; Marcondes et al. 2002, rat 4-day
  cycle, mouse values within ~10–20 % of these). The whole cycle is
  rotated so the proestrus midpoint anchors at 0.43, matching the
  human day-13 LH surge at `12 / 27 = 0.444`. Resulting positions:
  diestrus = 0.09, proestrus = 0.43, estrus = 0.55, metestrus = 0.71.
  Proestrus and estrus sit close together because the surge +
  ovulation window is only ~24 h of a ~4-day cycle; diestrus
  correspondingly occupies more than half the normalized axis.
  **Per-species rodent panels in §2.2 / §2.3 still use equal-width
  stages** (one tick per stage) for readability; the duration
  weighting is applied only to the cross-species panels here.
- **Spiny mouse:** explicit normalized phase positions from
  [cycle_seed.csv](data/seed/cycle_seed.csv) (menses 0.10,
  proliferative 0.40, secretory 0.75).

**Units are still per-species** so absolute magnitudes are not
directly comparable -- these plots are for **phase alignment**, not
amplitude comparison.

### 3.1 Estradiol

![Estradiol across species](plots/cycle_cross_species_estradiol.png)

The pre-ovulatory E2 peak lands at roughly the same normalized
position (~40–45 % through cycle) across human, mouse, rat, and
spiny mouse, supporting phase alignment but **not** evidence of
magnitude equivalence (human peak is 6–7× higher than mouse / rat
in absolute pg/mL; spiny mouse falls between human and rodent at
~116 pg/mL proliferative-phase mean). The spiny mouse is the only
non-human species in the seed with a documented **secretory-phase
E2 secondary elevation** (~107 pg/mL, lavender points in
Bellofiore 2021 Hum Reprod Fig 1e), giving it a human-like two-peak
shape rather than the single proestrus peak of mouse and rat (see
§2.4).

### 3.2 Progesterone

![Progesterone across species](plots/cycle_cross_species_progesterone.png)

This is the most biologically informative cross-species panel. Three
shapes are now visible:

- **Human** — broad luteal plateau (days 14–25, roughly 50–90 % of
  cycle), peak ~12 ng/mL.
- **Mouse / rat** — brief metestrus / post-ovulatory P4 bump only;
  no sustained luteal plateau.
- **Spiny mouse** — luteal-phase rise that sits **above** both the
  human plateau and rodent bump in absolute ng/mL (early follicular
  ~47 → late luteal ~135 ng/mL), but only two data points are
  seeded, so the *width* of the plateau is unresolved from this
  source. The 3-day menstrual / 6–10 day cycle implies the luteal
  phase spans roughly the last 25–35 % of the cycle, intermediate
  between the broad human plateau and the brief rodent bump.

The earlier hypothesis (luteal P4 *plateau width* as the candidate
correlate of spontaneous decidualization) is now joined by the
Bellofiore et al. 2018 alternative: the **luteal-to-follicular P4
fold-change** (≥3×) as the unique endocrine signature of
menstruating species. The spiny mouse 3.1× fold-change cleanly
separates muroid lineage from spontaneous decidualization —
**muroid mouse and rat are <3×; muroid spiny mouse is ≥3×** — so the
trait, not the clade, tracks with the fold-change. This is the one
result in the workspace that argues the spiny mouse is a high-value
comparator regardless of its phylogenetic position. Still
hypothesis-generating, not a tested claim.

### 3.3 LH and FSH

![LH across species](plots/cycle_cross_species_lh.png)
![FSH across species](plots/cycle_cross_species_fsh.png)

LH and FSH surges co-localize at the same normalized cycle position
(~40–45 %) in human, mouse, and rat. Shape is consistent; this is the
expected cross-species conservation and serves mainly as a sanity
check that the normalization is sensible. The spiny mouse is
**deliberately absent** from these two panels — see §2.4 for the
technical-ceiling reason (no validated assay + skin autotomy
precludes serial sampling). Adding spiny mouse traces here would
require new assay development.

### 3.4 Prolactin

![Prolactin across species](plots/cycle_cross_species_prolactin.png)

Two species are now on this panel: rat (ng/mL by RIA, proestrus
surge) and spiny mouse (IU/L by DXI immunoassay, late-luteal
secondary surge). Different assay families and different units — do
not compare amplitudes. The point of the panel is **phase
placement**: the rat prolactin surge co-localizes with the
pre-ovulatory LH surge, whereas the spiny mouse surge sits in the
late luteal phase (near 0.9 of normalized cycle), the same window
where the Bellofiore 2018 review reports a late-luteal prolactin
secondary surge in humans and rhesus. Adding a human prolactin
trace remains a clear next step before any further claim is made.

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
- **Cross-species rodent cycle normalization is now duration-weighted**
  using canonical rat 4-day-cycle hours (Marcondes et al. 2002):
  proestrus 12 h, estrus 12 h, metestrus 21 h, diestrus 57 h. This
  replaces the earlier equal-width approximation. Two residual
  caveats remain: (a) mouse durations are assumed equal to rat
  within ~10-20 %, not separately curated; (b) per-species rodent
  panels in §2.2 / §2.3 still draw stages at equal width for tick
  readability -- the duration weighting is applied only to the
  cross-species panels in §3.
- **The contraception matrix is endogenous-hormone-only.** Exogenous
  PK curves (ethinyl estradiol levels, levonorgestrel serum
  concentration vs time) are Phase 6B and deliberately deferred.
- **Primate rows are exploratory.** Do not let them carry
  quantitative claims. The spiny mouse row is plotted but seeded
  from only three assays (P4, prolactin, E2) with assay-family and
  unit caveats noted in §2.4 and §3.4; LH and FSH are unseedable
  for this species at present (see §2.4).
- **Spiny mouse P4 is at the source ceiling, not a curation gap.**
  Bellofiore 2018 (JME) Table 2 was inspected directly on
  2026-05-26 and confirmed to report P4 in only TWO phase bins
  (early follicular, late luteal), in contrast to the same paper's
  Table 3 prolactin which is reported in four bins. The companion
  Bellofiore 2021 (Hum Reprod) aging paper measures only estradiol
  and testosterone, not P4. The two-point spiny-mouse P4 trace
  therefore reflects what the published literature provides, not
  an incomplete seed.

## 6. Suggested next iterations (none committed)

In rough order of yield-per-effort:

1. **Species-specific rodent stage durations.** Current cross-species
   normalization uses rat 4-day-cycle hours (Marcondes 2002) for
   both mouse and rat. Splitting into separate mouse- vs rat-specific
   duration tables would tighten the per-species placement of estrus
   and metestrus by a few percent of the normalized axis. Probably
   not yield-positive unless a downstream claim depends on the
   distinction.
2. **Human prolactin series** to complete the cross-species prolactin
   panel.
3. **Contraception Phase 6B (exogenous PK curves)** -- only if
   regulatory or evolutionary questions explicitly need it; the
   endogenous-impact matrix is sufficient for atlas-context framing
   alone.

**Not on the next-iteration list:** spiny mouse LH and FSH absolute
serum values. As of 2021 these are not measurable in the species
(no validated immunoassay; skin autotomy precludes serial blood
draws). Closing that gap is a wet-lab methods project, not a
literature-fetch task, and is out of scope for this workspace.

None of these are commitments. Each would be a new scoped task,
proposed against PLAN.md, not a silent expansion of the year-one
atlas.
