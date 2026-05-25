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


def apply_species_overrides(
    set_name: str,
    mapped: list[str],
    species: str,
    kind: str,
) -> list[str]:
    """Augment a backbone-mapped gene list with per-species overrides.

    Parameters
    ----------
    set_name
        Gene-set key (e.g. ``"decidual_score"`` or ``"decidual_stromal"``).
    mapped
        Genes after ortholog-backbone mapping into ``species``.
    species
        Target species name (e.g. ``"mouse"``).
    kind
        One of ``"cell_type_markers"`` or ``"score_gene_sets"`` — matches the
        sub-key under ``species_overrides[species]`` in markers.yaml.

    Returns
    -------
    list[str]
        Augmented gene list with adds merged in and removes filtered out,
        preserving order and de-duplicating.
    """
    cfg = load_config("markers")
    overrides = cfg.get("species_overrides", {}).get(species, {}).get(kind, {}).get(set_name, {})
    if not overrides:
        return list(mapped)

    add = list(overrides.get("add") or [])
    remove = set(overrides.get("remove") or [])

    seen: set[str] = set()
    result: list[str] = []
    for g in list(mapped) + add:
        if g in remove or g in seen:
            continue
        seen.add(g)
        result.append(g)

    n_added = len(set(add) - set(mapped) - remove)
    n_removed = len(set(mapped) & remove)
    if n_added or n_removed:
        logger.info(
            "species_overrides[%s][%s][%s]: +%d added, -%d removed",
            species,
            kind,
            set_name,
            n_added,
            n_removed,
        )
    return result
