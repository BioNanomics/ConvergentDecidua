"""Normalize metadata fields to a standard vocabulary.

Provides mapping functions for species, assay, cycle stage, and
sample/donor metadata. Used by annotate.py to harmonize .obs columns.
"""

from __future__ import annotations

# -------------------------------------------------------------------------
# Species normalization
# -------------------------------------------------------------------------

_SPECIES_MAP = {
    "homo sapiens": "human",
    "human": "human",
    "hs": "human",
    "h. sapiens": "human",
    "mus musculus": "mouse",
    "mouse": "mouse",
    "mm": "mouse",
    "m. musculus": "mouse",
}


def normalize_species(raw: str) -> str:
    """Map a raw species string to the canonical name."""
    key = raw.strip().lower()
    if key in _SPECIES_MAP:
        return _SPECIES_MAP[key]
    return raw.strip().lower()


# -------------------------------------------------------------------------
# Assay normalization
# -------------------------------------------------------------------------

_ASSAY_MAP = {
    "scrna-seq": "scRNA-seq",
    "scrna": "scRNA-seq",
    "10x chromium": "scRNA-seq",
    "10x 3'": "scRNA-seq",
    "10x 5'": "scRNA-seq",
    "snrna-seq": "snRNA-seq",
    "snrna": "snRNA-seq",
    "scatac-seq": "scATAC-seq",
    "scatac": "scATAC-seq",
    "bulk-rna-seq": "bulk-RNA-seq",
    "bulk rna-seq": "bulk-RNA-seq",
    "rna-seq": "bulk-RNA-seq",
    "atac-seq": "ATAC-seq",
    "chip-seq": "ChIP-seq",
    "spatial": "spatial",
    "visium": "spatial",
}


def normalize_assay(raw: str) -> str:
    """Map a raw assay string to the canonical name."""
    key = raw.strip().lower()
    if key in _ASSAY_MAP:
        return _ASSAY_MAP[key]
    return raw.strip()


# -------------------------------------------------------------------------
# Cycle stage normalization
# -------------------------------------------------------------------------

_CYCLE_STAGE_MAP = {
    "proliferative": "proliferative",
    "early proliferative": "proliferative",
    "late proliferative": "proliferative",
    "early secretory": "early-secretory",
    "early-secretory": "early-secretory",
    "mid secretory": "mid-secretory",
    "mid-secretory": "mid-secretory",
    "late secretory": "late-secretory",
    "late-secretory": "late-secretory",
    "menstrual": "menstrual",
    "menses": "menstrual",
    "decidualized": "decidualized",
    "decidualization": "decidualized",
    "in vitro decidualized": "decidualized",
    "pregnancy": "pregnancy",
    "early pregnancy": "pregnancy",
    "implantation": "pregnancy",
    "diestrus": "diestrus",
    "estrus": "estrus",
    "proestrus": "proestrus",
    "metestrus": "metestrus",
}


def normalize_cycle_stage(raw: str) -> str:
    """Map a raw cycle stage string to the canonical name."""
    if not raw or raw.lower() in ("na", "nan", "unknown", ""):
        return "unknown"
    key = raw.strip().lower()
    if key in _CYCLE_STAGE_MAP:
        return _CYCLE_STAGE_MAP[key]
    return raw.strip().lower()


# -------------------------------------------------------------------------
# Donor/sample normalization
# -------------------------------------------------------------------------


def normalize_donor(raw: str, dataset: str) -> str:
    """Create a globally unique donor ID: {dataset}_{raw_id}."""
    clean = raw.strip().replace(" ", "_") if raw else "unknown"
    return f"{dataset}_{clean}"


def normalize_sample(raw: str, dataset: str) -> str:
    """Create a globally unique sample ID: {dataset}_{raw_id}."""
    clean = raw.strip().replace(" ", "_") if raw else "unknown"
    return f"{dataset}_{clean}"
