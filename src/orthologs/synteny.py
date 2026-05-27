"""Per-locus synteny / alignment check via the Ensembl REST API.

This complements ``src/orthologs/backbone.py`` (a static genome-wide
BioMart table) with a **per-locus alignment-quality** record for the
protected-core decidual panel. It answers the gate-C-style sanity
question: "for each anchor gene, does the partner species have a 1:1
ortholog with a documented pairwise alignment, and what is its
% identity / dN/dS?" — without anyone clicking through CGV.

Used by the ``wombat orthologs synteny-check`` CLI; output is written
to ``results/orthologs/synteny_at_core_loci.parquet`` and surfaced in
the ortholog report.
"""

from __future__ import annotations

import datetime as dt
import logging
import time
from pathlib import Path
from typing import Any

import pyarrow as pa
import pyarrow.parquet as pq

logger = logging.getLogger(__name__)

ENSEMBL_REST = "https://rest.ensembl.org"

# Map short species names (as used in configs/species.yaml and
# configs/markers.yaml::species_overrides) to Ensembl REST species
# identifiers. Extend here as the year-one scope widens.
_SPECIES_REST = {
    "human": "homo_sapiens",
    "mouse": "mus_musculus",
    "rat": "rattus_norvegicus",
    "spiny_mouse": "acomys_russatus",
}

# Output schema (pyarrow). Kept flat so it round-trips cleanly to
# parquet and pandas without nested type gymnastics.
_SCHEMA = pa.schema(
    [
        ("anchor_species", pa.string()),
        ("anchor_symbol", pa.string()),
        ("anchor_gene_id", pa.string()),
        ("target_species", pa.string()),
        ("target_symbol", pa.string()),
        ("target_gene_id", pa.string()),
        ("orthology_type", pa.string()),
        ("perc_id_source", pa.float64()),
        ("perc_id_target", pa.float64()),
        ("dn_ds", pa.float64()),
        ("alignment_present", pa.bool_()),
        ("fetched_at", pa.string()),
    ]
)


def _empty_row(
    anchor_species: str,
    anchor_symbol: str,
    target_species: str,
    *,
    note: str,
) -> dict[str, Any]:
    """Return a row representing 'no ortholog found' for this pair."""
    return {
        "anchor_species": anchor_species,
        "anchor_symbol": anchor_symbol,
        "anchor_gene_id": "",
        "target_species": target_species,
        "target_symbol": "",
        "target_gene_id": "",
        "orthology_type": note,
        "perc_id_source": float("nan"),
        "perc_id_target": float("nan"),
        "dn_ds": float("nan"),
        "alignment_present": False,
        "fetched_at": dt.datetime.now(dt.UTC).isoformat(timespec="seconds"),
    }


def _fetch_homology(
    anchor_species: str,
    symbol: str,
    target_species: str,
    *,
    timeout: float = 30.0,
) -> dict | None:
    """Query Ensembl REST for orthologues of one symbol in one target.

    Returns the parsed JSON dict, or ``None`` on a clean "no data"
    response. Raises ``requests.HTTPError`` on transport-level errors.
    """
    import requests

    rest_anchor = _SPECIES_REST.get(anchor_species, anchor_species)
    rest_target = _SPECIES_REST.get(target_species, target_species)

    url = f"{ENSEMBL_REST}/homology/symbol/{rest_anchor}/{symbol}"
    params = {
        "target_species": rest_target,
        "type": "orthologues",
        "sequence": "none",
    }
    # We omit `format=condensed` because it strips perc_id and dn_ds;
    # the default (full) payload keeps them and `sequence=none` already
    # drops the heavy alignment strings.
    headers = {"Accept": "application/json"}

    logger.debug("Ensembl REST: %s ?target=%s", url, rest_target)
    resp = requests.get(url, params=params, headers=headers, timeout=timeout)
    if resp.status_code == 400:
        # Ensembl returns 400 for "unknown symbol in this species" — treat
        # as a clean miss, not a transport failure.
        return None
    resp.raise_for_status()
    return resp.json()


def _lookup_symbol(
    gene_id: str,
    *,
    timeout: float = 30.0,
    _cache: dict[str, str] = {},  # noqa: B006 — intentional module-level memoization
) -> str:
    """Resolve an Ensembl gene ID to its display symbol.

    The ``/homology/symbol/`` endpoint omits ``target.symbol``, so we
    follow up with ``/lookup/id/{id}?expand=0`` to populate the symbol
    column. Results are memoized per-process to amortize cost across
    repeat invocations (e.g. multi-target runs).
    """
    if not gene_id:
        return ""
    if gene_id in _cache:
        return _cache[gene_id]

    import requests

    url = f"{ENSEMBL_REST}/lookup/id/{gene_id}"
    try:
        resp = requests.get(
            url,
            params={"expand": 0},
            headers={"Accept": "application/json"},
            timeout=timeout,
        )
        if resp.status_code != 200:
            _cache[gene_id] = ""
            return ""
        symbol = resp.json().get("display_name", "") or ""
    except Exception as exc:  # noqa: BLE001 — symbol lookup is best-effort
        logger.warning("Ensembl lookup failed for %s: %s", gene_id, exc)
        symbol = ""
    _cache[gene_id] = symbol
    return symbol


def _parse_homology(
    payload: dict | None,
    anchor_species: str,
    anchor_symbol: str,
    target_species: str,
) -> list[dict[str, Any]]:
    """Flatten one Ensembl REST homology payload into output rows."""
    fetched_at = dt.datetime.now(dt.UTC).isoformat(timespec="seconds")

    if not payload or not payload.get("data"):
        return [_empty_row(anchor_species, anchor_symbol, target_species, note="no_data")]

    rows: list[dict[str, Any]] = []
    for entry in payload["data"]:
        homologies = entry.get("homologies", [])
        if not homologies:
            rows.append(
                _empty_row(
                    anchor_species,
                    anchor_symbol,
                    target_species,
                    note="no_homology",
                )
            )
            continue
        for h in homologies:
            src = h.get("source", {}) or {}
            tgt = h.get("target", {}) or {}
            rows.append(
                {
                    "anchor_species": anchor_species,
                    "anchor_symbol": anchor_symbol,
                    "anchor_gene_id": src.get("id", "") or entry.get("id", ""),
                    "target_species": target_species,
                    "target_symbol": tgt.get("symbol", "") or "",
                    "target_gene_id": tgt.get("id", "") or "",
                    "orthology_type": h.get("type", ""),
                    "perc_id_source": float(src.get("perc_id", float("nan"))),
                    "perc_id_target": float(tgt.get("perc_id", float("nan"))),
                    "dn_ds": float(h.get("dn_ds") or float("nan")),
                    "alignment_present": h.get("type", "").startswith("ortholog"),
                    "fetched_at": fetched_at,
                }
            )
    return rows


def check_synteny(
    symbols: list[str],
    target_species: list[str],
    *,
    anchor_species: str = "human",
    rest_pause_s: float = 0.1,
) -> pa.Table:
    """Query Ensembl REST for per-locus orthology + alignment metadata.

    Parameters
    ----------
    symbols
        Anchor-species gene symbols to look up.
    target_species
        Short species names (e.g. ``["mouse", "rat"]``); mapped to
        Ensembl REST identifiers via ``_SPECIES_REST``.
    anchor_species
        Short anchor-species name. Defaults to ``"human"``.
    rest_pause_s
        Sleep between requests to stay within Ensembl REST rate limits
        (15 req/s soft cap; we go conservative at 10 req/s default).

    Returns
    -------
    pa.Table
        One row per (anchor_symbol, target_species, homology hit).
    """
    rows: list[dict[str, Any]] = []
    for symbol in symbols:
        for target in target_species:
            try:
                payload = _fetch_homology(anchor_species, symbol, target)
            except Exception as exc:  # noqa: BLE001 — surface, do not fail batch
                logger.warning("Ensembl REST failed for %s → %s: %s", symbol, target, exc)
                rows.append(
                    _empty_row(
                        anchor_species, symbol, target, note=f"rest_error:{type(exc).__name__}"
                    )
                )
                continue
            rows.extend(_parse_homology(payload, anchor_species, symbol, target))
            time.sleep(rest_pause_s)

    # Backfill target_symbol via /lookup/id/{id} — the homology endpoint
    # returns target.id but not target.symbol. Memoized inside
    # _lookup_symbol so repeat gene IDs cost one round-trip each.
    for r in rows:
        if r["target_gene_id"] and not r["target_symbol"]:
            r["target_symbol"] = _lookup_symbol(r["target_gene_id"])
            time.sleep(rest_pause_s)

    return pa.Table.from_pylist(rows, schema=_SCHEMA)


def load_protected_core() -> list[str]:
    """Load the protected_core gene list from configs/markers.yaml."""
    from wombat.config import load_config

    markers = load_config("markers")
    core = markers.get("protected_core") or []
    if not core:
        msg = "configs/markers.yaml has no protected_core list"
        raise ValueError(msg)
    return list(core)


def run_synteny_check(
    output_path: Path,
    *,
    symbols: list[str] | None = None,
    target_species: list[str] | None = None,
    anchor_species: str = "human",
) -> pa.Table:
    """End-to-end: fetch + write parquet. Returns the table for callers."""
    if symbols is None:
        symbols = load_protected_core()
    if target_species is None:
        target_species = ["mouse"]

    logger.info(
        "Synteny check: anchor=%s, targets=%s, %d symbols",
        anchor_species,
        target_species,
        len(symbols),
    )
    table = check_synteny(symbols, target_species, anchor_species=anchor_species)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, output_path)
    logger.info("Wrote %d rows → %s", table.num_rows, output_path)
    return table
