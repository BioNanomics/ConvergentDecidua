# ConvergentDecidua

A reproducible comparative atlas for the evolution of decidualization.

## Background

The evolution of spontaneous decidualization represents one of the most intriguing examples of convergent evolution in mammalian reproductive biology. Across multiple distantly related mammalian lineages — including higher primates, several bat species, elephant shrews, and a limited number of additional taxa — the endometrium evolved the capacity to undergo cyclical decidualization prior to embryo implantation. In these species, decidual transformation is initiated as part of the reproductive cycle itself rather than being induced directly by embryonic attachment signals, as occurs in the majority of mammals.

This transition fundamentally altered the temporal relationship between maternal tissues and embryonic invasion. Instead of responding reactively to implantation, the uterus began preparing proactively, establishing a hormonally regulated decidual state before the presence of an embryo was confirmed. In species exhibiting spontaneous decidualization, this cyclical preconditioning is tightly linked to menstruation, invasive placentation, and extensive maternal immune modulation. These traits are largely absent in closely related mammals with implantation-induced decidualization, suggesting that spontaneous decidualization evolved independently multiple times under similar selective pressures.

The repeated emergence of this phenotype raises a central evolutionary question: what genetic and regulatory changes shift decidualization from an embryo-triggered event to an internally timed cyclic program?

One prevailing hypothesis is that spontaneous decidualization evolved as a maternal adaptation to increasingly invasive trophoblast behavior. In this model, maternal tissues gained the ability to preemptively regulate implantation, constrain trophoblast invasion, and assess embryo quality before extensive placental integration occurred. Such a transition would require not merely changes in individual genes, but a rewiring of endocrine responsiveness, stromal cell differentiation programs, inflammatory signaling, and temporal regulatory networks governing the reproductive cycle.

The convergent appearance of spontaneous decidualization across phylogenetically distant mammals creates a powerful natural experiment for comparative genomics. If similar phenotypes evolved independently, then shared molecular signatures may reveal the core genetic architectures capable of generating cyclical decidual timing. These signatures may include changes in cis-regulatory elements, progesterone responsiveness, transcription factor binding networks, epigenetic regulation, noncoding RNAs, or alterations in developmental timing genes controlling endometrial stromal cell fate transitions.

This work investigates the genomic and regulatory basis of decidual timing by comparing species with spontaneous decidualization to closely related species retaining implantation-induced decidualization. Rather than focusing solely on genes associated with decidual identity, the emphasis here is on the evolution of timing itself: the transition from embryo-dependent activation to autonomous cyclical initiation. By identifying convergent regulatory changes across independently evolved menstruating lineages, it may be possible to uncover the minimal genetic circuitry required to transform decidualization into an anticipatory maternal program.

Ultimately, understanding how spontaneous decidualization evolved may illuminate broader principles governing evolutionary changes in developmental timing, maternal-fetal conflict, reproductive immunology, and the evolution of complex endocrine-regulated cellular states.

**CLI**: `wombat` · **Visualization**: DecidualAtlas (Streamlit) · **Current milestone**: MVR 0.1

## What this project does

ConvergentDecidua builds a cross-species single-cell atlas to investigate how spontaneous decidualization evolved independently in menstruating mammals. The pipeline ingests public scRNA-seq, scATAC-seq, and bulk RNA-seq datasets from human and mouse, maps orthologs, integrates stromal cell populations, and scores 8 decidualization-related gene modules.

For the full scientific background — hypothesis, species rationale, dataset targets, AI model strategy (DeciduaAI), and long-term roadmap — see [BACKGROUND.md](BACKGROUND.md).

## Current status (MVR 0.1)

The pipeline has been executed end-to-end on real data:

| Step | Result |
|---|---|
| **Ortholog backbone** | 25,439 human↔mouse pairs (16,168 Tier 1 + 9,271 Tier 2) via Ensembl Compara |
| **Datasets fetched** | GSE127918 (human scRNA), GSE111976 (human scRNA), GSE226429 (mouse bulk) |
| **QC** | GSE127918 → 9,292 cells; GSE111976 → 1,578 cells; GSE226429 → 6 samples |
| **Integration** | 9,065 human stromal cells via Harmony with UMAP embedding |
| **Cell types** | 4 subtypes: 8,466 fibroblast, 394 decidual, 129 pre-decidual, 76 senescent |
| **Scoring** | 8 modules scored; decidual_score highest in decidual_stromal (0.74) |
| **Reports** | Methods, coverage, QC, ortholog, scoring (heatmap + violin plots), manifest |

See [PLAN.md](PLAN.md) for detailed status, known gaps, and the implementation plan.

### Known gaps

- **No mouse scRNA integrated yet** — GSE226417 requires R/Seurat (RData format), E-MTAB-11491 has 644 individual files
- **No scATAC data** — GSE183771 is a 7GB+ tar archive, deferred
- **Integration is human-only** — blocked by the mouse data gap above

## Repository structure

```
wombat/              CLI and orchestration (Click commands, config loader)
src/                 Analysis modules
  ingest/              GEO/ArrayExpress download → AnnData
  metadata/            .obs harmonization
  qc/                  scRNA, scATAC, bulk QC pipelines
  orthologs/           Ensembl BioMart / Compara ortholog mapping
  cell_states/         Marker-based annotation + Harmony/scVI integration
  scoring/             8 decidualization modules via scanpy.tl.score_genes
  reports/             Methods, coverage, QC, manifest generation
configs/             YAML configs (datasets, species, markers)
decidual_atlas/      Streamlit visualization app
workflows/           Snakemake rules
tests/               pytest suite
results/             Pipeline outputs (.gitignored)
  orthologs/           backbone.parquet, orthogroups
  raw/                 Downloaded source files
  processed/           Standardized h5ad per dataset
  qc/                  QC-filtered h5ad
  integrated/          Harmony/scVI joint embeddings
  scored/              Decidualization-scored h5ad
  reports/             Generated reports and figures
```

## Installation

```bash
# Clone
git clone https://github.com/BioNanomics/ConvergentDecidua.git
cd ConvergentDecidua

# Create virtual environment
python -m venv .venv
source .venv/bin/activate

# Install with all optional dependencies
pip install -e ".[all]"

# Or install specific groups
pip install -e ".[dev]"       # pytest, ruff, pre-commit
pip install -e ".[atlas]"     # streamlit, plotly
pip install -e ".[ingest]"    # GEOparse
pip install -e ".[workflow]"  # snakemake
```

Requires Python ≥3.9, <3.13. Additional runtime dependencies: `tabulate`, `harmonypy`.

## Wombat CLI

`wombat` is the command-line interface that drives the entire pipeline. Each step can be run independently or chained via Snakemake.

```
wombat [OPTIONS] COMMAND [ARGS]

Options:
  -v, --verbose    Increase verbosity (-v, -vv)

Commands:
  init              Validate that all required configs exist
  validate-config   Validate all YAML configuration files
  build-registry    Export dataset registry to Parquet and CSV
  fetch             Download datasets and convert to standardized AnnData
  qc                Run quality control on processed datasets
  orthologs build   Build ortholog backbone and orthogroup tables
  integrate         Integrate datasets across species
  score-decidua     Compute decidualization scores across all 8 modules
  generate-reports  Generate all pipeline reports and release manifest
  serve-atlas       Launch the DecidualAtlas Streamlit app
```

### Typical workflow

```bash
# 1. Validate configuration
wombat validate-config

# 2. Build ortholog backbone (human↔mouse)
wombat orthologs build

# 3. Fetch and QC datasets
wombat fetch --dataset GSE127918
wombat fetch --all-datasets
wombat qc --species human
wombat qc --species mouse

# 4. Integrate stromal cells
wombat integrate --mode stromal --method harmony

# 5. Score decidualization modules
wombat score-decidua

# 6. Generate reports
wombat generate-reports

# 7. Launch interactive viewer
wombat serve-atlas --port 8501
```

## Decidualization scoring modules

The pipeline scores each cell on 8 gene-set modules defined in `configs/markers.yaml`:

| Module | What it captures |
|---|---|
| `decidual_score` | Core decidualization signature (PRL, IGFBP1, FOXO1, etc.) |
| `progesterone_response_score` | Progesterone receptor pathway activity |
| `estrogen_response_score` | Estrogen receptor pathway activity |
| `stress_response_score` | Oxidative/cellular stress response |
| `senescence_score` | Cellular senescence markers |
| `immune_interface_score` | Stromal–immune crosstalk genes |
| `ECM_remodeling_score` | Extracellular matrix remodeling |
| `angiogenesis_score` | Angiogenesis and vascular remodeling |

## Configuration

All pipeline parameters are defined in YAML files under `configs/`:

- **`datasets.yaml`** — Dataset registry (accession, species, assay, priority)
- **`species.yaml`** — Species metadata (genome build, Ensembl release, decidualization mode)
- **`markers.yaml`** — Cell-type markers and decidualization gene sets

Config loading uses `wombat.config.load_config(name)` — never hardcode paths.

## How we got here

The project was built in two phases:

### Phase 1: Code scaffolding (Phases 1–10)

Starting from an empty repo, 10 implementation phases built the full pipeline code:

1. **Project skeleton** — pyproject.toml, Click CLI, CI, pre-commit
2. **Configuration** — YAML schemas for datasets, species, markers
3. **Data ingestion** — GEO + ArrayExpress downloaders → standardized AnnData
4. **Metadata harmonization** — Consistent `.obs` columns across all datasets
5. **QC pipelines** — scRNA (doublets, mito%, HVG), scATAC (TSS, LSI), bulk
6. **Ortholog mapping** — BioMart + Compara FTP with cross-validation
7. **Cell-state integration** — Marker scoring, Harmony/scVI joint embedding
8. **Decidualization scoring** — 8 gene-set modules via scanpy
9. **DecidualAtlas** — Streamlit app with UMAP, gene explorer, species comparison
10. **Reports** — Methods, coverage, QC, ortholog, manifest with checksums

### Phase 2: Execution against real data (E1–E5)

When the code hit real APIs and data, several fixes were needed:

- **BioMart mirrors down** → added Ensembl Compara FTP fallback for ortholog mapping
- **GEO URL prefix bug** → fixed accession slicing for 6+ digit accessions
- **DGE format detection** → broadened file pattern matching for tab-separated count matrices
- **HVG numpy incompatibility** → added `cell_ranger` fallback when `seurat_v3` fails
- **Sparse matrix h5ad write** → convert to dense for small bulk datasets
- **Harmony shape mismatch** → auto-detect `Z_corr` transposition
- **Python 3.9 compat** → `datetime.UTC` → `timezone.utc`, HVG flavor fallback

Each fix was committed individually and the pipeline re-run to validate.

## Development

```bash
# Lint
ruff check .

# Test
pytest

# Validate configs
wombat validate-config
```

CI runs lint, test, and config validation on every push via GitHub Actions.

## License

MIT
