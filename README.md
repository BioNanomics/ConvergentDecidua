# ConvergentDecidua

> A reproducible comparative atlas for the evolution of decidualization.

| Component | Name |
|---|---|
| **Repo** | `ConvergentDecidua` |
| **CLI** | `wombat` |
| **AI layer** | `DeciduaAI` (Phase 2) |
| **Database** | `DeciduaForge` (Phase 3) |
| **Visualization** | `DecidualAtlas` (Phase 4) |

**Current milestone**: [MVR 0.1](PLAN.md) — human + mouse atlas foundation.

Core public anchors include human cycle scRNA/scATAC datasets like GSE127918 and GSE183771, mouse decidualization datasets, and AI tools such as scGPT, Geneformer, LINGER, and Enformer. GSE127918 maps human stromal decidual trajectories by single-cell RNA-seq, while GSE183771 profiles human endometrial chromatin accessibility across the menstrual cycle. ([NCBI][1])

---

# ConvergentDecidua

## Mission

Build a reproducible comparative decidualization atlas to investigate the evolutionary origins of spontaneous decidualization, identify conserved gene networks underlying an emergent biological oscillator that helps control cycle length and decidualization, and detect convergent regulatory changes across humans, bats, spiny mouse, and ordinary laboratory mouse as an outgroup.

The project should produce:

1. a curated comparative dataset registry,
2. reproducible processed single-cell / bulk / chromatin objects,
3. ortholog and synteny mappings,
4. decidualization state labels,
5. regulatory network models,
6. sequence-level candidate regulatory elements,
7. ranked convergence hypotheses,
8. a database and visualization layer.

---

# 1. Scientific objective

## Primary question

What regulatory changes plausibly enabled spontaneous decidualization to evolve independently in menstruating mammals?

## Focal species

Use this initial logic:

| Species     | Role                                                                    |
| ----------- | ----------------------------------------------------------------------- |
| Human       | richest spontaneous decidualization reference                           |
| Bat         | independent spontaneous decidualization / menstruation lineage          |
| Spiny mouse | menstruating rodent; sequence-first unless more uterine omics are found |
| Lab mouse   | non-spontaneous decidualization outgroup                                |

## Working hypothesis

Spontaneous decidualization is not caused by one “switch gene,” but by regulatory rewiring of stromal-cell hormone, stress-response, immune, and extracellular-matrix programs.

---

# 2. Project phases

## Phase 1 — Atlas foundation

Goal: build a clean, reproducible comparative atlas.

Major workstreams:

```text
1. Dataset registry
2. Data ingestion
3. Metadata harmonization
4. Ortholog mapping
5. Cell-type harmonization
6. Decidualization scoring
7. Human/mouse/bat/spiny baseline comparison
```

## Phase 2 — DeciduaAI

Goal: use models to reduce search space.

Models to include:

| Model              | Use                                                   |
| ------------------ | ----------------------------------------------------- |
| scVI / scANVI      | batch correction, integration, label transfer         |
| scGPT              | cell-state embeddings and transfer learning           |
| Geneformer         | human-centric regulator prioritization                |
| GENIE3 / dynGENIE3 | baseline GRN inference                                |
| LINGER             | GRN inference when RNA+ATAC/multiome support exists   |
| Enformer / Borzoi  | sequence-to-regulation scoring                        |
| AlphaGenome        | optional frontier rescoring if access/licensing works |

scGPT is a single-cell foundation model trained on over 33 million cells and supports annotation, integration, perturbation-response prediction, and gene network inference. ([Nature][2]) Geneformer is another relevant foundation model for context-aware single-cell transcriptomics. ([Nature][3]) LINGER is useful for inferring gene regulatory networks from paired single-cell RNA/chromatin data. ([Nature][4]) Enformer predicts regulatory activity from long DNA sequence context and is useful for cis-regulatory evolution questions. ([Nature][5])

## Phase 3 — DeciduaForge

Goal: persist all processed evidence in a queryable database.

Store:

```text
datasets
samples
species
genes
orthologs
cell types
cell states
marker scores
GRN edges
candidate regulatory elements
sequence-model scores
convergence scores
```

Recommended database:

```text
PostgreSQL + DuckDB/Parquet + object storage
```

PostgreSQL for metadata and queries; Parquet/Zarr/H5AD for matrices.

## Phase 4 — DecidualAtlas

Goal: interactive exploration.

Views:

```text
species comparison
cell-state browser
decidualization marker heatmaps
gene/regulon explorer
candidate element explorer
trajectory viewer
convergence scoreboard
```

Stack:

```text
Streamlit or Dash
Plotly
DuckDB
AnnData/Zarr
```

---

# 3. Initial dataset targets

## Human

| Accession   | Type       | Use                                              |
| ----------- | ---------- | ------------------------------------------------ |
| GSE111976   | scRNA-seq  | human endometrium across natural menstrual cycle |
| GSE127918   | scRNA-seq  | decidual pathway / stromal trajectory            |
| GSE183771   | scATAC-seq | chromatin across menstrual cycle                 |
| GSE205477   | RNA-seq    | TRIM28 perturbation in human stromal cells       |
| GSE205473   | ATAC-seq   | TRIM28 perturbation accessibility                |
| GSE62475    | ChIP-seq   | progesterone receptor cistrome                   |
| E-MTAB-9260 | spatial    | spatial validation layer                         |

## Mouse

| Accession    | Type         | Use                                  |
| ------------ | ------------ | ------------------------------------ |
| E-MTAB-11491 | scRNA-seq    | cycling and decidualizing mouse FRT  |
| GSE226417    | scRNA-seq    | early pregnancy decidua / uterus     |
| GSE226429    | bulk RNA-seq | in vitro decidualization time course |
| GSE122376    | bulk RNA-seq | mouse decidualization models         |
| GSE205480    | scRNA-seq    | perturbation reference               |
| E-MTAB-12105 | spatial      | mouse spatial validation             |

## Bat

| Resource        | Type         | Use                                 |
| --------------- | ------------ | ----------------------------------- |
| GCF_021234435.1 | genome       | Jamaican fruit bat reference        |
| PRJNA1251670    | snRNA-seq    | maternal-fetal interface validation |
| PRJNA1251235    | single-cell  | optional bat validation             |
| PRJNA1251203    | bulk RNA-seq | pseudobulk validation               |

## Spiny mouse

| Resource        | Type          | Use                                  |
| --------------- | ------------- | ------------------------------------ |
| GCA_029890205.1 | genome        | sequence-first regulatory comparison |
| PRJNA182705     | transcriptome | annotation / ortholog support        |

The spiny mouse genome assembly GCA_029890205.1 is available as a chromosome-level draft assembly. ([National Genomics Data Center][6])

---

# 4. Repo structure

```text
ConvergentDecidua/
  README.md
  pyproject.toml
  environment.yml
  docker/
    Dockerfile
  configs/
    datasets.yaml
    species.yaml
    markers.yaml
    models.yaml
    database.yaml
  workflows/
    Snakefile
    rules/
      fetch.smk
      qc.smk
      orthologs.smk
      integrate.smk
      grn.smk
      sequence.smk
      convergence.smk
      reports.smk
  wombat/
    cli.py
    config.py
    logging.py
  src/
    ingest/
    metadata/
    qc/
    orthologs/
    cell_states/
    scoring/
    grn/
    sequence/
    convergence/
    reports/
  decidua_ai/
    scvi_runner.py
    scgpt_runner.py
    geneformer_runner.py
    genie3_runner.py
    linger_runner.py
    enformer_runner.py
  decidua_forge/
    schema.sql
    models.py
    loaders.py
    queries.py
  decidual_atlas/
    app.py
    pages/
  notebooks/
  tests/
  docs/
  results/
```

---

# 5. Wombat CLI commands

```bash
wombat init
wombat validate-config
wombat fetch --dataset GSE127918
wombat build-registry
wombat qc --species human
wombat orthologs build
wombat integrate --mode stromal
wombat score-decidua
wombat infer-grn --method genie3
wombat infer-grn --method linger
wombat score-sequence --model enformer
wombat find-convergence
wombat load-forge
wombat serve-atlas
```

---

# 6. Core data model

## Dataset registry

```yaml
- accession: GSE127918
  species: human
  tissue: endometrium
  assay: scRNA-seq
  condition: decidualization_timecourse
  spontaneous_decidualization: true
  menstruates: true
  role: discovery
  priority: high
```

## Gene table

```text
gene_id
symbol
species
ensembl_id
orthogroup_id
human_anchor_gene
orthology_type
```

## Cell-state table

```text
cell_id
dataset
species
sample
donor
cycle_stage
cell_type
cell_state
decidual_score
progesterone_score
stress_score
senescence_score
```

## Candidate regulatory element table

```text
element_id
species
chrom
start
end
nearest_gene
linked_gene
motifs
te_family
accessibility_score
sequence_model_score
convergence_score
```

---

# 7. Ortholog strategy

Use three tiers:

## Tier 1 — strict one-to-one orthologs

Use for primary cross-species integration.

## Tier 2 — orthogroups

Use for one-to-many and many-to-many gene families.

## Tier 3 — regulatory synteny

Use for enhancers and noncoding regions where basewise alignment may fail.

Primary tools:

```text
Ensembl Compara / BioMart
g:Profiler / g:Orth
OrthoDB
BLAST / minimap2 fallback
UCSC chain/liftOver where available
```

---

# 8. Cell-state harmonization

Initial ontology:

```text
stromal_fibroblast
pre_decidual_stromal
decidual_stromal
senescent_decidual
perivascular
epithelial_glandular
epithelial_luminal
endothelial
immune_uNK
macrophage
T_cell
other
```

Initial marker set:

```text
PGR
ESR1
FOXO1
HOXA10
HAND2
GATA2
PRL
IGFBP1
IL15
WNT4
BMP2
LEFTY2
MMPs
CXCLs
```

---

# 9. Decidualization scoring

Build score modules:

```text
decidual_score
progesterone_response_score
estrogen_response_score
stress_response_score
senescence_score
immune_interface_score
ECM_remodeling_score
angiogenesis_score
```

Each score should be:

```text
gene_set
species_mapped_gene_set
per-cell score
per-sample pseudobulk score
confidence score
```

---

# 10. DeciduaAI: model-assisted narrowing

## Objective

Reduce the search space from thousands of genes/elements to a ranked list of candidate convergent regulators.

## Input

```text
integrated stromal AnnData
ortholog table
human/mouse/bat/spiny genome sequences
ATAC peaks where available
GRN priors
decidualization scores
```

## Output

```text
candidate_regulators.parquet
candidate_elements.bed
candidate_network_edges.parquet
convergence_scores.parquet
```

## Model workflow

1. **scVI/scANVI**

   * integrate datasets
   * label transfer
   * species/donor-aware latent space

2. **scGPT**

   * learn stromal-state embeddings
   * test whether bat/mouse/spiny-like states project toward human decidual trajectories

3. **Geneformer**

   * run human stromal cells through perturbation/ranking workflows
   * prioritize likely upstream regulators

4. **GENIE3 / dynGENIE3**

   * infer expression-only GRNs
   * use as baseline and sanity check

5. **LINGER**

   * use where chromatin + expression support exists
   * infer TF → regulatory element → target gene networks

6. **Enformer / Borzoi**

   * score candidate regulatory regions
   * run motif ablations
   * run sequence swaps between human/bat/spiny/mouse

7. **Convergence ranking**

   * combine expression, GRN, chromatin, sequence, synteny, motif, and TE evidence.

---

# 11. Convergence scoring

Candidate score should combine:

```text
ortholog confidence
decidual-state specificity
human evidence
mouse outgroup contrast
bat recurrence
spiny mouse sequence support
GRN centrality
chromatin support
motif gain/loss
TE-family enrichment
sequence-model effect
literature support
```

Example output:

```text
rank
gene
regulator_type
human_score
bat_score
spiny_score
mouse_outgroup_score
grn_score
sequence_score
convergence_score
recommended_validation
```

---

# 12. DeciduaForge database

## Purpose

Make the atlas queryable and reproducible.

## Suggested stack

```text
PostgreSQL
SQLAlchemy
DuckDB
Parquet
Zarr
AnnData
Alembic migrations
```

## Key tables

```text
species
datasets
samples
genes
orthologs
cell_types
cell_states
gene_sets
scores
grn_edges
regulatory_elements
sequence_scores
convergence_candidates
workflow_runs
```

## Example query

```sql
SELECT gene_symbol, convergence_score, grn_score, sequence_score
FROM convergence_candidates
WHERE cell_state = 'decidual_stromal'
ORDER BY convergence_score DESC
LIMIT 100;
```

---

# 13. DecidualAtlas visualization layer

## Views

1. Dataset registry browser
2. Species comparison dashboard
3. Cell-state UMAP / embedding viewer
4. Decidualization score heatmap
5. Gene explorer
6. GRN explorer
7. Regulatory element browser
8. Convergence leaderboard
9. Evidence card per candidate

## Candidate evidence card

For each gene/regulatory element:

```text
gene / element
species evidence
ortholog status
cell-state specificity
expression pattern
GRN edges
chromatin support
motif changes
sequence-model prediction
TE annotation
validation priority
```

---

# 14. Implementation plan

See [PLAN.md](PLAN.md) for the phased implementation plan, dependency graph, MVR 0.1 dataset selections, and per-epic step breakdowns.

### Epic summary

| Epic | Scope | MVR |
|---|---|---|
| A | Project skeleton (repo, pyproject, Docker, Snakemake, CLI, CI) | 0.1 |
| B | Dataset registry | 0.1 |
| C | Data ingestion | 0.1 |
| D | Metadata harmonization | 0.1 |
| E | QC | 0.1 |
| F | Orthologs and synteny | 0.1 |
| G | Cell-state harmonization | 0.1 |
| H | Scoring | 0.1 |
| I | DeciduaAI model runners | 0.3 |
| J | Convergence engine | 0.3 |
| K | DeciduaForge database | 1.0 |
| L | DecidualAtlas visualization | 0.1 (baseline), 1.0 (full) |
| M | Reports and reproducibility | 0.1 (baseline), 1.0 (full) |

---

# 15. Minimum viable releases

| Release | Scope |
|---|---|
| **MVR 0.1** | Human + mouse, processed matrices, ortholog backbone, stromal cell-state harmonization, decidualization scoring, baseline atlas |
| **MVR 0.2** | Add bat, spiny mouse genome, sequence windows, GRN baseline |
| **MVR 0.3** | Add DeciduaAI (scGPT/Geneformer/GENIE3/Enformer), ranked convergence candidates |
| **MVR 1.0** | DeciduaForge database, full DecidualAtlas, reproducible workflow, candidate evidence cards, paper-style report |

---

# 16. Definition of success

A successful first version should produce this statement:

```text
Given public human, mouse, bat, and spiny mouse resources, ConvergentDecidua identifies a ranked set of stromal regulators and candidate regulatory elements whose expression, network position, chromatin support, sequence-model scores, and outgroup contrast are consistent with convergent evolution of spontaneous decidualization.
```

It does **not** need to prove causality. It should generate high-quality, reproducible hypotheses for follow-up CRISPRi, reporter, organoid, or stromal-cell perturbation experiments.

[1]: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE127918&utm_source=chatgpt.com "GEO Accession viewer"
[2]: https://www.nature.com/articles/s41592-024-02201-0?utm_source=chatgpt.com "scGPT: toward building a foundation model for single-cell multi-omics using generative AI | Nature Methods"
[3]: https://www.nature.com/articles/s41592-024-02305-7?utm_source=chatgpt.com "Large-scale foundation model on single-cell transcriptomics | Nature Methods"
[4]: https://www.nature.com/articles/s41587-024-02182-7?utm_source=chatgpt.com "Inferring gene regulatory networks from single-cell multiome data using atlas-scale external data | Nature Biotechnology"
[5]: https://www.nature.com/articles/s41592-021-01252-x?utm_source=chatgpt.com "Effective gene expression prediction from sequence by integrating long-range interactions | Nature Methods"
[6]: https://ngdc.cncb.ac.cn/gwh/ncbi_assembly/36959/show?utm_source=chatgpt.com "Genome Warehouse"
