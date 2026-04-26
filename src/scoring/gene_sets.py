"""Load and manage gene sets from configs/markers.yaml.

Implements all 8 decidualization scoring modules (README §9):
  1. decidual_score
  2. progesterone_response_score
  3. estrogen_response_score
  4. stress_response_score
  5. senescence_score
  6. immune_interface_score
  7. ECM_remodeling_score
  8. angiogenesis_score
"""

from __future__ import annotations

import logging

from wombat.config import load_config

logger = logging.getLogger(__name__)


def load_score_gene_sets() -> dict[str, list[str]]:
    """Load all scoring gene sets from markers.yaml.

    Returns
    -------
    dict[str, list[str]]
        Mapping of score name → list of human gene symbols.
    """
    cfg = load_config("markers")
    gene_sets = cfg.get("score_gene_sets", {})

    if not gene_sets:
        msg = "No score_gene_sets found in configs/markers.yaml"
        raise ValueError(msg)

    logger.info(
        "Loaded %d scoring modules: %s",
        len(gene_sets),
        ", ".join(gene_sets.keys()),
    )
    return gene_sets


def load_cell_type_markers() -> dict[str, list[str]]:
    """Load cell-type marker gene sets from markers.yaml.

    Returns
    -------
    dict[str, list[str]]
        Mapping of cell type → list of human marker genes.
    """
    cfg = load_config("markers")
    markers = cfg.get("cell_type_markers", {})
    return {k: v for k, v in markers.items() if v}
