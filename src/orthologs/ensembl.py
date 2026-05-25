"""Fetch orthologs from Ensembl BioMart REST API.

Queries the Ensembl BioMart for one-to-one orthologs between two species.
Results are cached locally as Parquet to avoid redundant API calls.
"""

from __future__ import annotations

import logging
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

# BioMart dataset names per species
_DATASETS = {
    "human": "hsapiens_gene_ensembl",
    "mouse": "mmusculus_gene_ensembl",
}

# BioMart attribute names for ortholog queries (source → target)
_ORTHOLOG_ATTRS = {
    "mouse": {
        "homolog_ensembl_gene": "mmusculus_homolog_ensembl_gene",
        "homolog_associated_gene_name": "mmusculus_homolog_associated_gene_name",
        "homolog_orthology_type": "mmusculus_homolog_orthology_type",
        "homolog_orthology_confidence": "mmusculus_homolog_orthology_confidence",
    },
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
        Source species name (must be in _DATASETS).
    target : str
        Target species name (must be in _ORTHOLOG_ATTRS).
    cache_dir : Path, optional
        Directory to cache raw results. If a cache file exists, it is loaded
        instead of querying the API.

    Returns
    -------
    pa.Table
        Table with columns: human_gene_id, human_symbol, mouse_gene_id,
        mouse_symbol, orthology_type, confidence.
    """
    import requests

    if source not in _DATASETS:
        msg = f"Unknown source species: {source}. Available: {list(_DATASETS)}"
        raise ValueError(msg)
    if target not in _ORTHOLOG_ATTRS:
        msg = f"Unknown target species: {target}. Available: {list(_ORTHOLOG_ATTRS)}"
        raise ValueError(msg)

    # Check cache
    if cache_dir is not None:
        cache_path = cache_dir / f"ensembl_{source}_{target}_raw.parquet"
        if cache_path.exists():
            logger.info("Loading cached Ensembl orthologs from %s", cache_path)
            return pq.read_table(cache_path)

    # Build and send query — try mirrors on failure
    xml = _build_query_xml(_DATASETS[source], _ORTHOLOG_ATTRS[target])

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
            resp = r
            break  # success
        except requests.RequestException as exc:
            logger.warning("Mirror %s failed: %s", mirror_url, exc)
            continue

    if resp is None:
        logger.warning("All BioMart mirrors failed. Falling back to Ensembl Compara FTP...")
        table = _fetch_compara_ftp(source, target)
        logger.info("Retrieved %d ortholog rows from Compara FTP", len(table))

        if cache_dir is not None:
            cache_dir.mkdir(parents=True, exist_ok=True)
            pq.write_table(table, cache_path)
            logger.info("Cached to %s", cache_path)

        return table

    # Parse TSV response
    table = _parse_biomart_response(resp.text, source, target)
    logger.info("Retrieved %d ortholog rows from Ensembl", len(table))

    # Cache result
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        pq.write_table(table, cache_path)
        logger.info("Cached to %s", cache_path)

    return table


def _parse_biomart_response(text: str, source: str, target: str) -> pa.Table:
    """Parse BioMart TSV response into a PyArrow Table."""
    import csv

    reader = csv.DictReader(StringIO(text), delimiter="\t")
    rows = list(reader)

    if not rows:
        msg = "BioMart returned no data"
        raise RuntimeError(msg)

    # Standardize column names
    source_gene_id_col = "Gene stable ID"
    source_symbol_col = "Gene name"

    # The target columns use BioMart naming — find them dynamically
    fieldnames = reader.fieldnames or []

    # Map BioMart response columns to our standard names
    # BioMart returns e.g. "Mouse gene stable ID", "Mouse gene name",
    # "Mouse homology type", "Mouse orthology confidence [0 low, 1 high]"
    target_gene_id_col = _find_column(fieldnames, "mouse", "gene", "stable")
    target_symbol_col = _find_column(fieldnames, "mouse", "gene", "name")
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

_COMPARA_SPECIES = {
    "human": "homo_sapiens",
    "mouse": "mus_musculus",
}


def _fetch_compara_ftp(source: str, target: str) -> pa.Table:
    """Download pre-built ortholog table from Ensembl Compara FTP.

    The Compara FTP provides TSV files with all homology types. We filter
    for the target species and return the same schema as BioMart.
    """
    import gzip
    import io

    import requests

    src_species = _COMPARA_SPECIES.get(source)
    tgt_species = _COMPARA_SPECIES.get(target)
    if not src_species or not tgt_species:
        msg = f"Unknown species for Compara FTP: {source} → {target}"
        raise ValueError(msg)

    # Download from target species directory (contains target→source mappings)
    url = f"{_COMPARA_FTP}/{tgt_species}/"
    logger.info("Listing Compara FTP directory: %s", url)
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()

    # Find the protein homology file
    import re

    matches = re.findall(r'href="(Compara\.\d+\.protein_default\.homologies\.tsv\.gz)"', resp.text)
    if not matches:
        msg = f"No Compara protein homology file found at {url}"
        raise RuntimeError(msg)

    file_url = f"{_COMPARA_FTP}/{tgt_species}/{matches[0]}"
    logger.info("Downloading Compara orthologs from %s", file_url)
    resp = requests.get(file_url, timeout=300, stream=True)
    resp.raise_for_status()

    # Decompress and parse — filter for source species
    raw = gzip.decompress(resp.content)
    reader = io.StringIO(raw.decode("utf-8"))

    # TSV columns: gene_stable_id, protein_stable_id, species, identity,
    #   homology_type, homology_gene_stable_id, homology_protein_stable_id,
    #   homology_species, homology_identity, dn, ds, goc_score,
    #   wga_coverage, is_high_confidence, homology_id
    import csv

    tsv_reader = csv.DictReader(reader, delimiter="\t")

    source_ids = []
    source_symbols = []  # Compara doesn't have symbols — filled later
    target_ids = []
    target_symbols = []
    orth_types = []
    confidences = []

    for row in tsv_reader:
        if row.get("homology_species") != src_species:
            continue

        target_ids.append(row["gene_stable_id"])
        target_symbols.append("")  # Compara TSV doesn't include gene symbols
        source_ids.append(row["homology_gene_stable_id"])
        source_symbols.append("")
        orth_types.append(row.get("homology_type", ""))
        conf = row.get("is_high_confidence", "0")
        confidences.append(int(conf) if conf.isdigit() else 0)

    logger.info("Parsed %d %s→%s ortholog rows from Compara", len(source_ids), source, target)

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
