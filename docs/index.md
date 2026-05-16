# ConvergentDecidua Documentation

Reproducible comparative single-cell atlas for the evolution of decidualization across mammals.

**Repository**: [github.com/BioNanomics/ConvergentDecidua](https://github.com/BioNanomics/ConvergentDecidua)

---

## Overview

ConvergentDecidua builds a cross-species single-cell atlas to investigate how spontaneous decidualization evolved independently in menstruating mammals. The pipeline ingests public scRNA-seq, scATAC-seq, and bulk RNA-seq datasets from human and mouse, maps orthologs, integrates stromal cell populations, and scores decidualization-related gene modules.

## Key Components

- **`wombat` CLI** — Command-line interface for all pipeline steps (Click-based)
- **DecidualAtlas** — Streamlit visualization app for exploring the integrated atlas
- **Snakemake workflows** — Reproducible data processing rules
- **Scoring engine** — Generic gene-set scoring framework for 8 decidualization modules

## Current Milestone: MVR 0.1

| Step | Result |
|---|---|
| Ortholog backbone | 25,439 human↔mouse pairs via Ensembl Compara |
| Datasets | GSE127918, GSE111976 (human scRNA), GSE226429 (mouse bulk) |
| Integration | Harmony on 9,065 human stromal cells, 4 subtypes identified |
| Scoring | 8 modules scored (decidual, progesterone response, senescence, etc.) |

## Links

- [Background & Scientific Context](https://github.com/BioNanomics/ConvergentDecidua/blob/main/BACKGROUND.md)
- [Implementation Plan](https://github.com/BioNanomics/ConvergentDecidua/blob/main/PLAN.md)
- [Reports](https://github.com/BioNanomics/ConvergentDecidua/tree/main/results/reports)

## Citation

If you use this software, please cite using the [CITATION.cff](https://github.com/BioNanomics/ConvergentDecidua/blob/main/CITATION.cff) file.

## License

[MIT](https://github.com/BioNanomics/ConvergentDecidua/blob/main/LICENSE)
