# Comparative species matrix

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
  `exploratory` species are present in the matrix and source manifest
  but only plotted if a downloadable source is confirmed;
  `reference_only` species inform the comparative framing but are not
  plotted.

| species | closeness | spont. decidu. | menstruation | genome | hormone data | downloadability | first-pass |
|---|---|---|---|---|---|---|---|
| human | reference | yes | yes | T2T complete | excellent, daily serum | structured table | plotted |
| rhesus macaque | close Old World monkey | yes | yes | T2T complete | good, cycle resolved | manual extraction | reference only |
| common marmoset | close New World monkey | unconfirmed, likely no | no | chromosome scale | moderate, ovarian cycle | manual extraction | exploratory |
| gray mouse lemur | strepsirrhine | likely no | no | chromosome scale | sparse | reference only | exploratory |
| spiny mouse (*Acomys*) | rodent (Acomys) | yes | yes | chromosome scale | sparse, phase resolved | manual extraction | plotted |
| mouse | rodent outgroup | no | no | T2T complete | good, stage resolved | manual extraction | plotted |
| rat | rodent outgroup | no | no | chromosome scale | good, stage resolved | manual extraction | plotted |

## Notes

- Rhesus macaque is intentionally **reference only** here because it
  sits on the wrong side of the spontaneous-decidualization axis for
  the trait-negative control role, even though it has strong cycle
  hormone data. It is the right comparator for menstruation work but
  the wrong comparator for "closely related, without spontaneous
  decidualization."
- Marmoset is marked **unconfirmed, likely no**: its endometrial
  literature describes cycle-dependent uterine priming around
  implantation without a clean cyclical spontaneous decidual program.
  This is the residual question from the prior research turn.
- Mouse lemur is the cleaner candidate for a phenotype-clean primate
  trait-negative control on the spontaneous-decidualization axis, but
  endocrine data is too sparse for first-pass plotting.
- Spiny mouse is now **plotted** as a positive control on the
  menstruation/spontaneous-decidualization axis using Bellofiore et
  al. 2017 (AJOG) / 2018 (JME) Tables 2-3 for P4 and prolactin, plus
  Bellofiore et al. 2021 (Human Reproduction) Fig 1e for estradiol
  (~50 pg/mL menses nadir, ~140 pg/mL proliferative peak, secondary
  elevation in the secretory phase). **LH and FSH absolute serum
  values cannot be measured in this species at present** -- no
  validated cross-reactive immunoassay exists and skin autotomy
  prevents the serial blood draws that would be needed for cycle-
  resolved gonadotropin tracking (McKenna 2020 PLoS ONE; McKenna
  2021 Sci Rep). This is recorded as a species-level technical
  ceiling in [sources.yaml](sources.yaml)
  (`spiny_mouse_gonadotropin_ceiling`), not a missing-paper gap.
- Bat is intentionally absent from this matrix because the project
  scope downgraded it to sequence-first; revisit if a downloadable
  bat hormone-cycle table appears.

See [sources.yaml](sources.yaml) for source-level provenance behind
each `hormone data` rating.
