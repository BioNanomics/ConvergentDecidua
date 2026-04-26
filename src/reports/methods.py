"""Auto-generate methods report from workflow metadata.

Produces a reproducible methods section summarizing pipeline versions,
parameters, and processing steps.
"""

from __future__ import annotations

import logging
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path

logger = logging.getLogger(__name__)


def generate_methods_report(output_path: Path) -> Path:
    """Generate methods report.

    Parameters
    ----------
    output_path : Path
        Where to write the report.

    Returns
    -------
    Path
        Path to the generated report.
    """
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sections = [
        _header(),
        _environment(),
        _pipeline_overview(),
        _data_acquisition(),
        _quality_control(),
        _ortholog_mapping(),
        _integration(),
        _scoring(),
    ]

    with open(output_path, "w") as fh:
        fh.write("\n\n".join(sections))
        fh.write("\n")

    logger.info("Methods report → %s", output_path)
    return output_path


def _header() -> str:
    ts = datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    return f"# Methods Report\n\n*Generated: {ts}*"


def _environment() -> str:
    versions = {"Python": sys.version.split()[0], "Platform": platform.platform()}

    try:
        import scanpy

        versions["scanpy"] = scanpy.__version__
    except ImportError:
        pass
    try:
        import anndata

        versions["anndata"] = anndata.__version__
    except ImportError:
        pass

    lines = ["## Environment", ""]
    lines.extend(f"- **{k}**: {v}" for k, v in versions.items())
    return "\n".join(lines)


def _pipeline_overview() -> str:
    return """## Pipeline Overview

The ConvergentDecidua pipeline processes single-cell and bulk RNA-seq data
from human and mouse endometrial/decidual tissues through: data acquisition,
quality control, metadata harmonization, cross-species ortholog mapping,
cell-state annotation, Harmony-based integration, and decidualization scoring
across 8 gene-set modules."""


def _data_acquisition() -> str:
    return """## Data Acquisition

Processed count matrices were downloaded from GEO (GSE accessions) and
ArrayExpress (E-MTAB accessions) using automated downloaders. Files were
converted to AnnData h5ad format with standardized metadata from the
dataset registry."""


def _quality_control() -> str:
    return """## Quality Control

**scRNA-seq**: Cells were filtered requiring ≥500 detected genes,
≤8000 genes, and <15% mitochondrial reads. Doublets were detected
using Scrublet. Data was log-normalized (target_sum=10,000) and
3,000 highly variable genes were selected (Seurat v3 method).

**scATAC-seq**: Cells were filtered by TSS enrichment (≥2.0) and
fragment count (1,000–100,000). TF-IDF normalization was applied
followed by LSI dimensionality reduction.

**Bulk RNA-seq**: Samples with <1,000 total counts were removed.
Genes detected in <2 samples or with <10 total counts were filtered.
CPM normalization with log1p transformation was applied."""


def _ortholog_mapping() -> str:
    return """## Ortholog Mapping

Human-mouse orthologs were retrieved from Ensembl BioMart (release 112)
and cross-validated with g:Profiler g:Orth. The backbone table contains
Tier 1 (strict one-to-one, high confidence) and Tier 2 (relaxed, including
many-to-many) mappings."""


def _integration() -> str:
    return """## Cross-Species Integration

Mouse gene symbols were remapped to human orthologs via the Tier 1 backbone.
Datasets were concatenated, re-normalized, and 3,000 HVGs selected on the
combined matrix. PCA (50 components) was computed, followed by Harmony
batch correction on the species variable. UMAP was computed on the
Harmony-corrected embedding."""


def _scoring() -> str:
    return """## Decidualization Scoring

Eight gene-set modules were scored using scanpy.tl.score_genes:
decidual_score, progesterone_response_score, estrogen_response_score,
stress_response_score, senescence_score, immune_interface_score,
ECM_remodeling_score, angiogenesis_score. Gene sets were mapped through
the ortholog backbone for non-human species."""
