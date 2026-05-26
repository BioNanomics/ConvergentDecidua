# Workflow DAG

The Snakemake DAG for the cross-species stromal atlas pipeline.

## Pipeline overview

```mermaid
flowchart TD
    Datasets[(configs/datasets.yaml)]
    Fetch[fetch_dataset<br/>wombat fetch]
    QC[qc_dataset<br/>wombat qc]
    Orthologs[build_orthologs<br/>wombat orthologs build]
    Integrate[integrate_stromal<br/>wombat integrate --mode stromal]
    Reports[generate_reports<br/>wombat generate-reports]

    Datasets --> Fetch
    Fetch --> QC
    QC --> Integrate
    Orthologs --> Integrate
    Integrate --> Reports

    classDef io fill:#eef,stroke:#447,color:#000
    classDef step fill:#efe,stroke:#474,color:#000
    class Datasets io
    class Fetch,QC,Orthologs,Integrate,Reports step
```

## Per-rule artifacts

| Rule | Output | Producer |
|---|---|---|
| `fetch_dataset` | `results/processed/{accession}.h5ad` | `wombat fetch` |
| `qc_dataset` | `results/qc/{accession}.h5ad` | `wombat qc` |
| `build_orthologs` | `results/orthologs/backbone.parquet` | `wombat orthologs build` |
| `integrate_stromal` | `results/integrated/stromal_cross_species.h5ad` | `wombat integrate` |
| `generate_reports` | `results/reports/*.md` | `wombat generate-reports` |

## Reproducing the DAG

The full DOT-format DAG is checked in at [dag.dot](dag.dot) and can be
re-generated with:

```bash
snakemake --snakefile workflows/Snakefile --dag --forceall > docs/dag.dot
# Optional render (requires graphviz):
dot -Tsvg docs/dag.dot > docs/dag.svg
```

CI dry-runs the DAG on every push (`validate-workflow` job) to catch
broken `include:` paths and rule syntax errors without running data
pipelines.
