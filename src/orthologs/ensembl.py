"""Fetch orthologs from Ensembl BioMart REST API.

Queries the Ensembl BioMart for one-to-one orthologs between two species.
Results are cached locally as Parquet to avoid redundant API calls.

Dataset names, BioMart attribute prefixes, and Compara FTP species
identifiers are resolved per-call from ``configs/species.yaml`` via the
``ensembl_dataset``, ``ensembl_prefix``, and ``ensembl_species`` fields,
so adding a new species is a config edit only.
"""

from __future__ import annotations

import logging
from functools import lru_cache
from io import StringIO
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

logger = logging.getLogger(__name__)

BIOMART_URL = "https://www.ensembl.org/biomart/martservice"

# Fallback mirrors in case primary is down
_BIOMART_MIRRORS = [
    "https://www.ensembl.org/biomart/martservice",
    "https://useast.ensembl.org/biomart/martservice",
    "https://asia.ensembl.org/biomart/martservice",
]


@lru_cache(maxsize=1)
def _species_index() -> dict[str, dict]:
    """Return ``{species_name: entry_dict}`` from configs/species.yaml."""
    from wombat.config import load_config

    entries = load_config("species")
    return {entry["name"]: entry for entry in entries}


def _species_field(name: str, field: str) -> str:
    """Look up a required field on a species entry; raise informatively."""
    idx = _species_index()
    if name not in idx:
        msg = f"Unknown species: {name!r}. Known: {sorted(idx)}"
        raise ValueError(msg)
    entry = idx[name]
    if field not in entry or entry[field] is None:
        msg = (
            f"species.yaml entry for {name!r} is missing required field "
            f"{field!r} (needed for Ensembl ortholog lookup)."
        )
        raise ValueError(msg)
    return entry[field]


def _ortholog_attrs(target_prefix: str) -> dict[str, str]:
    """Build the four BioMart ortholog attribute names for a target species."""
    return {
        "homolog_ensembl_gene": f"{target_prefix}_homolog_ensembl_gene",
        "homolog_associated_gene_name": f"{target_prefix}_homolog_associated_gene_name",
        "homolog_orthology_type": f"{target_prefix}_homolog_orthology_type",
        "homolog_orthology_confidence": f"{target_prefix}_homolog_orthology_confidence",
    }


def _build_query_xml(
    source_dataset: str,
    target_attrs: dict[str, str],
) -> str:
    """Build BioMart XML query for ortholog retrieval."""
    attr_lines = "\n".join(
        [
            '      <Attribute name="ensembl_gene_id" />',
            '      <Attribute name="external_gene_name" />',
        ]
        + [f'      <Attribute name="{v}" />' for v in target_attrs.values()]
    )

    return f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1">
  <Dataset name="{source_dataset}" interface="default">
{attr_lines}
  </Dataset>
</Query>"""


def fetch_ensembl_orthologs(
    source: str = "human",
    target: str = "mouse",
    cache_dir: Path | None = None,
) -> pa.Table:
    """Fetch one-to-one orthologs from Ensembl BioMart.

    Parameters
    ----------
    source : str
        Source species name (must exist in ``configs/species.yaml``).
    target : str
        Target species name (must exist in ``configs/species.yaml`` with
        an ``ensembl_prefix`` field).
    cache_dir : Path, optional
        Directory to cache raw results. If a cache file exists, it is loaded
        instead of querying the API.

    Returns
    -------
    pa.Table
        Table with columns: ``{source}_gene_id``, ``{source}_symbol``,
        ``{target}_gene_id``, ``{target}_symbol``, ``orthology_type``,
        ``confidence``.
    """
    import requests

    source_dataset = _species_field(source, "ensembl_dataset")
    target_prefix = _species_field(target, "ensembl_prefix")
    target_attrs = _ortholog_attrs(target_prefix)

    # Check cache
    if cache_dir is not None:
        cache_path = cache_dir / f"ensembl_{source}_{target}_raw.parquet"
        if cache_path.exists():
            logger.info("Loading cached Ensembl orthologs from %s", cache_path)
            return pq.read_table(cache_path)

    # Build and send query — try mirrors on failure
    xml = _build_query_xml(source_dataset, target_attrs)

    resp = None
    for mirror_url in _BIOMART_MIRRORS:
        logger.info("Querying BioMart (%s) for %s→%s orthologs...", mirror_url, source, target)
        try:
            r = requests.get(mirror_url, params={"query": xml}, timeout=300)
            r.raise_for_status()
            if r.text.startswith("Query ERROR"):
                logger.warning(
                    "Mirror %s returned BioMart query error: %s",
                    mirror_url,
                    r.text[:500],
                )
                continue
            # Detect HTML responses (e.g. Ensembl "Service unavailable" page,
            # 403 Forbidden bodies, captive portals). A real TSV response
            # starts with the column header "Gene stable ID\t...".
            ctype = (r.headers.get("content-type") or "").lower()
            stripped = r.text.lstrip()
            if "html" in ctype or stripped.startswith(("<", "<!DOCTYPE", "<!doctype")):
                logger.warning(
                    "Mirror %s returned an HTML page (likely maintenance); first 200 chars: %s",
                    mirror_url,
                    stripped[:200].replace("\n", " "),
                )
                continue
            resp = r
            break  # success
        except requests.RequestException as exc:
            logger.warning("Mirror %s failed: %s", mirror_url, exc)
            continue

    if resp is None:
        logger.warning("All BioMart mirrors failed. Falling back to Ensembl Compara FTP...")
        table = _fetch_compara_ftp(source, target)
        logger.info("Retrieved %d ortholog rows from Compara FTP", len(table))

        if len(table) == 0:
            msg = (
                f"All BioMart mirrors failed AND Compara FTP returned 0 rows "
                f"for {source}->{target}. Refusing to cache an empty result; "
                f"retry when Ensembl is healthy."
            )
            raise RuntimeError(msg)

        if cache_dir is not None:
            cache_dir.mkdir(parents=True, exist_ok=True)
            pq.write_table(table, cache_path)
            logger.info("Cached to %s", cache_path)

        return table

    # Parse TSV response
    table = _parse_biomart_response(resp.text, source, target)
    logger.info("Retrieved %d ortholog rows from Ensembl", len(table))

    if len(table) == 0:
        # The BioMart response was valid TSV (passed the HTML / Query ERROR
        # guards) but contained zero usable ortholog rows. This usually means
        # the upstream stream was truncated mid-response. Don't cache the
        # empty parse — let the next run retry.
        msg = (
            f"BioMart parsed 0 ortholog rows for {source}->{target}. "
            f"Likely a truncated response; refusing to cache. Retry."
        )
        raise RuntimeError(msg)

    # Cache result
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        pq.write_table(table, cache_path)
        logger.info("Cached to %s", cache_path)

    return table


def _parse_biomart_response(text: str, source: str, target: str) -> pa.Table:
    """Parse BioMart TSV response into a PyArrow Table.

    BioMart returns ``Gene stable ID`` / ``Gene name`` for the source
    species (no species qualifier) and ``<Target> gene stable ID`` /
    ``<Target> gene name`` / ``<Target> homology type`` /
    ``<Target> orthology confidence`` for the target. The exact target
    label depends on Ensembl's display name for the species; we therefore
    pick the target columns as "the stable-ID/name column that is NOT the
    bare source one" so the parser works for any target species without
    knowing its display label.
    """
    import csv

    reader = csv.DictReader(StringIO(text), delimiter="\t")
    rows = list(reader)

    if not rows:
        msg = "BioMart returned no data"
        raise RuntimeError(msg)

    fieldnames = reader.fieldnames or []

    source_gene_id_col = "Gene stable ID"
    source_symbol_col = "Gene name"
    target_gene_id_col = _find_other_column(fieldnames, "stable", exclude=source_gene_id_col)
    target_symbol_col = _find_other_column(
        fieldnames, "name", exclude=source_symbol_col, also_exclude=(source_gene_id_col,)
    )
    orthology_type_col = _find_column(fieldnames, "homology", "type")
    confidence_col = _find_column(fieldnames, "confidence")

    source_ids = []
    source_symbols = []
    target_ids = []
    target_symbols = []
    orth_types = []
    confidences = []

    for row in rows:
        src_id = row.get(source_gene_id_col, "").strip()
        tgt_id = row.get(target_gene_id_col, "").strip()
        if not src_id or not tgt_id:
            continue

        source_ids.append(src_id)
        source_symbols.append(row.get(source_symbol_col, "").strip())
        target_ids.append(tgt_id)
        target_symbols.append(row.get(target_symbol_col, "").strip())
        orth_types.append(row.get(orthology_type_col, "").strip())
        confidences.append(int(row.get(confidence_col, "0").strip() or "0"))

    return pa.table(
        {
            f"{source}_gene_id": source_ids,
            f"{source}_symbol": source_symbols,
            f"{target}_gene_id": target_ids,
            f"{target}_symbol": target_symbols,
            "orthology_type": orth_types,
            "confidence": confidences,
        }
    )


def _find_other_column(
    fieldnames: list[str],
    keyword: str,
    *,
    exclude: str,
    also_exclude: tuple[str, ...] = (),
) -> str:
    """Find the column containing ``keyword`` that is not ``exclude``.

    Used to locate the target-species ``... stable ID`` / ``... name``
    columns without needing to know the species' display label.
    """
    skip = {exclude, *also_exclude}
    for col in fieldnames:
        if col in skip:
            continue
        if keyword.lower() in col.lower():
            return col
    msg = f"Could not find a column containing {keyword!r} other than {exclude!r} in {fieldnames}"
    raise KeyError(msg)


def _find_column(fieldnames: list[str], *keywords: str) -> str:
    """Find a column name containing all keywords (case-insensitive)."""
    for col in fieldnames:
        col_lower = col.lower()
        if all(kw.lower() in col_lower for kw in keywords):
            return col
    msg = f"Could not find column matching keywords: {keywords} in {fieldnames}"
    raise KeyError(msg)


# Ensembl Compara FTP — pre-built ortholog tables (fallback when BioMart is down)
_COMPARA_FTP = "https://ftp.ensembl.org/pub/current_tsv/ensembl-compara/homologies"


def _fetch_compara_ftp(source: str, target: str) -> pa.Table:
    """Download pre-built ortholog table from Ensembl Compara FTP.

    The Compara FTP provides TSV files with all homology types. We filter
    for the target species and return the same schema as BioMart. Latin
    binomial identifiers come from the ``ensembl_species`` field in
    ``configs/species.yaml``.
    """
    import gzip
    import io

    import requests

    src_species = _species_field(source, "ensembl_species")
    tgt_species = _species_field(target, "ensembl_species")

    import csv
    import re

    def _try_dir(dir_species: str, want_other: str) -> pa.Table:
        """Fetch the Compara homologies TSV under ``dir_species`` and filter
        rows where ``homology_species == want_other``.

        Ensembl Compara only ships each pair on one side (typically the
        smaller / less-paired genome's directory). For some species
        (e.g. ``papio_anubis``) the per-species file does not include
        ``homo_sapiens`` at all, so we may need to try both directions.
        """
        url = f"{_COMPARA_FTP}/{dir_species}/"
        logger.info("Listing Compara FTP directory: %s", url)
        resp = requests.get(url, timeout=60)
        resp.raise_for_status()
        matches = re.findall(
            r'href="(Compara\.\d+\.protein_default\.homologies\.tsv\.gz)"', resp.text
        )
        if not matches:
            msg = f"No Compara protein homology file found at {url}"
            raise RuntimeError(msg)
        file_url = f"{url}{matches[0]}"
        logger.info("Downloading Compara orthologs from %s", file_url)
        resp = requests.get(file_url, timeout=600)
        resp.raise_for_status()
        raw = gzip.decompress(resp.content)
        reader = io.StringIO(raw.decode("utf-8"))
        tsv_reader = csv.DictReader(reader, delimiter="\t")

        # Identify source vs target columns by which side carries which species.
        # The file under dir_species has species=dir_species, homology_species=other.
        # We want source on the source side regardless of directory.
        source_ids: list[str] = []
        target_ids: list[str] = []
        orth_types: list[str] = []
        confidences: list[int] = []

        if dir_species == src_species:
            src_col, tgt_col = "gene_stable_id", "homology_gene_stable_id"
        else:
            src_col, tgt_col = "homology_gene_stable_id", "gene_stable_id"

        for row in tsv_reader:
            if row.get("homology_species") != want_other:
                continue
            source_ids.append(row[src_col])
            target_ids.append(row[tgt_col])
            orth_types.append(row.get("homology_type", ""))
            conf = row.get("is_high_confidence", "0")
            confidences.append(int(conf) if conf.isdigit() else 0)

        return pa.table(
            {
                f"{source}_gene_id": source_ids,
                f"{source}_symbol": [""] * len(source_ids),
                f"{target}_gene_id": target_ids,
                f"{target}_symbol": [""] * len(target_ids),
                "orthology_type": orth_types,
                "confidence": confidences,
            }
        )

    # Try target directory first (usually smaller), fall back to source directory.
    table = _try_dir(tgt_species, src_species)
    if len(table) == 0:
        logger.warning(
            "Compara file under %s/ contained no %s homologies; retrying via %s/ directory.",
            tgt_species,
            src_species,
            src_species,
        )
        table = _try_dir(src_species, tgt_species)

    logger.info("Parsed %d %s→%s ortholog rows from Compara", len(table), source, target)

    return table
