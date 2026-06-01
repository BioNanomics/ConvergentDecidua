"""Per-species transcript-to-gene mapping from Ensembl BioMart.

The GSE155170 deposit (and a handful of other older bulk RNA-seq
deposits) ships transcript-level count matrices with ``ENS<spp>T...``
identifiers in ``var_names``. The Tier-B ortholog backbones are
indexed on **gene** IDs (``ENS<spp>G...``), so a transcript→gene
rollup is required before any ortholog mapping can happen.

This module:

1. Queries Ensembl BioMart for the species's ``ensembl_transcript_id``,
   ``ensembl_gene_id``, and ``external_gene_name`` attributes, caching
   the result to ``results/orthologs/cache/tx2gene_<species>.parquet``.
2. Aggregates an AnnData whose ``var_names`` are transcript IDs into a
   gene-level AnnData by summing counts of all transcripts that map to
   the same gene.

BioMart access reuses the 3-mirror retry + HTML-maintenance guard from
``src/orthologs/ensembl.py`` (kept inline here rather than imported to
avoid coupling the ingest layer to the orthologs module's internals).
"""

from __future__ import annotations

import logging
from io import StringIO
from pathlib import Path

import anndata as ad
import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
import scipy.sparse as sp

logger = logging.getLogger(__name__)

_BIOMART_MIRRORS = (
    "https://www.ensembl.org/biomart/martservice",
    "https://useast.ensembl.org/biomart/martservice",
    "https://asia.ensembl.org/biomart/martservice",
)


def _build_tx2gene_xml(dataset: str, chromosomes: str | None = None) -> str:
    """Build the BioMart query XML, optionally restricted to ``chromosomes``.

    ``chromosomes`` is a comma-separated list of ``chromosome_name``
    values; when provided, the query is filtered to those regions so a
    large transcriptome can be pulled in small chunks that stay under
    server-side proxy timeouts.
    """
    chrom_filter = (
        f'      <Filter name="chromosome_name" value="{chromosomes}" />\n' if chromosomes else ""
    )
    return f"""<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1">
  <Dataset name="{dataset}" interface="default">
{chrom_filter}      <Attribute name="ensembl_transcript_id" />
      <Attribute name="ensembl_gene_id" />
      <Attribute name="external_gene_name" />
  </Dataset>
</Query>"""


def _parse_tx2gene_tsv(text: str) -> pa.Table:
    import csv

    reader = csv.DictReader(StringIO(text), delimiter="\t")
    fields = reader.fieldnames or []
    # BioMart's column labels for these attributes:
    #   "Transcript stable ID", "Gene stable ID", "Gene name"
    tx_col = next((f for f in fields if "Transcript" in f and "stable ID" in f), None)
    gene_col = next((f for f in fields if f == "Gene stable ID"), None)
    sym_col = next((f for f in fields if f == "Gene name"), None)
    if tx_col is None or gene_col is None:
        msg = f"Unexpected BioMart columns for tx2gene: {fields!r}"
        raise RuntimeError(msg)

    tx_ids: list[str] = []
    gene_ids: list[str] = []
    symbols: list[str] = []
    for row in reader:
        t = (row.get(tx_col) or "").strip()
        g = (row.get(gene_col) or "").strip()
        if not t or not g:
            continue
        tx_ids.append(t)
        gene_ids.append(g)
        symbols.append((row.get(sym_col) or "").strip() if sym_col else "")
    return pa.table(
        {
            "transcript_id": tx_ids,
            "gene_id": gene_ids,
            "gene_symbol": symbols,
        }
    )


# Connect fast so a dead mirror is abandoned quickly; allow a longer
# read budget once the connection is live and data is streaming.
_BIOMART_CONNECT_TIMEOUT = 15
_BIOMART_READ_TIMEOUT = 180
# Emit a heartbeat at most this often while streaming, so a slow but
# alive transfer is visibly making progress instead of looking hung.
_BIOMART_HEARTBEAT_SECONDS = 5.0


def _fetch_biomart_tsv(
    xml: str,
    species: str,
    mirrors: tuple[str, ...] = _BIOMART_MIRRORS,
) -> tuple[str, list[str]]:
    """Download a BioMart TSV, streaming with progress + flaky reporting.

    Tries each mirror in turn. A mirror is considered *flaky* (but not
    fatal) when it times out, errors, returns an HTML maintenance page,
    or returns a BioMart ``Query ERROR``; the issue is recorded and the
    next mirror is tried. While a healthy mirror streams data, a
    heartbeat is logged at most every ``_BIOMART_HEARTBEAT_SECONDS`` so
    a slow-but-alive transfer never looks hung.

    ``mirrors`` overrides the default mirror list; pass a single-host
    tuple (e.g. an Ensembl archive ``.../biomart/martservice``) to pin
    the query to an annotation release that matches an older deposit.

    Returns
    -------
    (text, flaky) : the decoded TSV body and a list of human-readable
    notes describing any mirrors that misbehaved before success.

    Raises
    ------
    RuntimeError if every mirror fails.
    """
    import time

    import requests

    flaky: list[str] = []
    for mirror_url in mirrors:
        start = time.monotonic()
        logger.info("Querying BioMart (%s) for %s tx2gene...", mirror_url, species)
        try:
            with requests.get(
                mirror_url,
                params={"query": xml},
                timeout=(_BIOMART_CONNECT_TIMEOUT, _BIOMART_READ_TIMEOUT),
                stream=True,
            ) as r:
                r.raise_for_status()
                ctype = (r.headers.get("content-type") or "").lower()
                chunks: list[bytes] = []
                total = 0
                last_beat = start
                for chunk in r.iter_content(chunk_size=64 * 1024):
                    if not chunk:
                        continue
                    chunks.append(chunk)
                    total += len(chunk)
                    now = time.monotonic()
                    if now - last_beat >= _BIOMART_HEARTBEAT_SECONDS:
                        logger.info(
                            "  ...%s: %.1f MB received (%.0fs elapsed)",
                            species,
                            total / 1e6,
                            now - start,
                        )
                        last_beat = now
                text = b"".join(chunks).decode("utf-8", errors="replace")
            elapsed = time.monotonic() - start
        except requests.RequestException as exc:
            note = f"{mirror_url} failed after {time.monotonic() - start:.0f}s: {exc}"
            logger.warning(note)
            flaky.append(note)
            continue

        stripped = text.lstrip()
        if text.startswith("Query ERROR"):
            note = f"{mirror_url} returned BioMart error: {text[:200]}"
            logger.warning(note)
            flaky.append(note)
            continue
        if "html" in ctype or stripped.startswith(("<", "<!DOCTYPE", "<!doctype")):
            note = f"{mirror_url} returned an HTML maintenance page"
            logger.warning("%s; trying next mirror.", note)
            flaky.append(note)
            continue

        logger.info(
            "BioMart %s OK for %s: %.1f MB in %.0fs",
            mirror_url,
            species,
            total / 1e6,
            elapsed,
        )
        return text, flaky

    msg = (
        f"All BioMart mirrors failed for tx2gene({species}). "
        f"Compara FTP does not host transcript→gene mappings; "
        f"retry when Ensembl is healthy. Details: {'; '.join(flaky)}"
    )
    raise RuntimeError(msg)


_ENSEMBL_REST = "https://rest.ensembl.org"
# Main karyotype regions are fetched one chunk at a time (they hold the
# bulk of the genes); small unplaced scaffolds are batched together so
# the request count stays modest.
_SCAFFOLD_BATCH = 50
# A chunk can fail when every mirror briefly blinks out (e.g. the primary
# flips to a maintenance page for a few seconds). Retry the chunk with
# backoff rather than aborting the whole multi-chunk run.
_CHUNK_RETRY_BACKOFF = (5, 15, 30, 60)


def _fetch_chromosome_names(ensembl_species: str) -> list[str]:
    """List top-level region names (chromosomes + scaffolds) via Ensembl REST.

    Using the assembly's full top-level region set — not just the named
    karyotype — guarantees chunked fetches still cover genes on unplaced
    scaffolds, so the concatenated result matches a single-shot query.
    """
    import requests

    url = f"{_ENSEMBL_REST}/info/assembly/{ensembl_species}"
    r = requests.get(
        url,
        params={"content-type": "application/json"},
        timeout=(_BIOMART_CONNECT_TIMEOUT, 60),
    )
    r.raise_for_status()
    data = r.json()
    regions = [str(reg["name"]) for reg in data.get("top_level_region", []) if reg.get("name")]
    if not regions:
        msg = f"Ensembl REST returned no top-level regions for {ensembl_species}."
        raise RuntimeError(msg)
    return regions


def _group_chromosome_chunks(regions: list[str]) -> list[str]:
    """Group region names into BioMart ``chromosome_name`` filter values.

    Karyotype-style regions (short names like ``1``..``19``, ``X``,
    ``Y``, ``MT``) are emitted one per chunk; longer scaffold names are
    grouped ``_SCAFFOLD_BATCH`` at a time into comma-separated values.
    """
    main: list[str] = []
    scaffolds: list[str] = []
    for name in regions:
        if len(name) <= 2 or name.upper() in {"MT", "X", "Y"}:
            main.append(name)
        else:
            scaffolds.append(name)
    chunks: list[str] = [m for m in main]
    for i in range(0, len(scaffolds), _SCAFFOLD_BATCH):
        chunks.append(",".join(scaffolds[i : i + _SCAFFOLD_BATCH]))
    return chunks


def _fetch_chunk_with_retry(
    xml: str,
    label: str,
    mirrors: tuple[str, ...] = _BIOMART_MIRRORS,
) -> tuple[str, list[str]]:
    """Fetch one chunk, retrying with backoff if all mirrors blink out.

    Transient outages (the primary mirror flipping to a maintenance page
    for a few seconds, a one-off 504) should not abort a long multi-chunk
    run, so a chunk that exhausts every mirror is retried after a short
    wait before giving up.
    """
    import time

    flaky: list[str] = []
    attempts = len(_CHUNK_RETRY_BACKOFF) + 1
    for attempt in range(1, attempts + 1):
        try:
            text, chunk_flaky = _fetch_biomart_tsv(xml, label, mirrors=mirrors)
            flaky.extend(chunk_flaky)
            return text, flaky
        except RuntimeError as exc:
            flaky.append(str(exc))
            if attempt == attempts:
                raise
            wait = _CHUNK_RETRY_BACKOFF[attempt - 1]
            logger.warning(
                "Chunk %s failed (attempt %d/%d); retrying in %ds...",
                label,
                attempt,
                attempts,
                wait,
            )
            time.sleep(wait)
    # Unreachable: the loop either returns or raises.
    raise RuntimeError(f"Exhausted retries for chunk {label}")  # pragma: no cover


def _fetch_tx2gene_chunked(
    dataset: str,
    ensembl_species: str,
    species: str,
    mirrors: tuple[str, ...] = _BIOMART_MIRRORS,
) -> tuple[pa.Table, list[str]]:
    """Fetch a tx2gene table in per-chromosome chunks and concatenate.

    Each chunk is a small query that completes well under server-side
    proxy timeouts, sidestepping the 504 that a single full-transcriptome
    query triggers for large genomes (e.g. mouse).
    """
    import pandas as pd

    regions = _fetch_chromosome_names(ensembl_species)
    chunks = _group_chromosome_chunks(regions)
    logger.info(
        "Chunked tx2gene for %s: %d regions → %d query chunks",
        species,
        len(regions),
        len(chunks),
    )

    frames: list[pd.DataFrame] = []
    flaky: list[str] = []
    for i, chrom in enumerate(chunks, start=1):
        label = (
            chrom if len(chrom) <= 20 else f"{chrom[:17]}... ({chrom.count(',') + 1} scaffolds)"
        )
        logger.info("  chunk %d/%d (%s) for %s", i, len(chunks), label, species)
        xml = _build_tx2gene_xml(dataset, chromosomes=chrom)
        text, chunk_flaky = _fetch_chunk_with_retry(xml, f"{species}:{label}", mirrors=mirrors)
        flaky.extend(chunk_flaky)
        frames.append(_parse_tx2gene_tsv(text).to_pandas())

    combined = pd.concat(frames, ignore_index=True)
    combined = combined.drop_duplicates("transcript_id", ignore_index=True)
    table = pa.Table.from_pandas(combined, preserve_index=False)
    return table, flaky


def _optional_species_field(species: str, field: str) -> str | None:
    """Return an optional species-config field, or None if it is absent."""
    from src.orthologs.ensembl import _species_field

    try:
        return _species_field(species, field)
    except ValueError:
        return None


def fetch_tx2gene(
    species: str,
    cache_dir: Path | None = None,
    chunked: bool = False,
    biomart_host: str | None = None,
) -> pa.Table:
    """Fetch a transcript→gene table for ``species`` from Ensembl BioMart.

    Parameters
    ----------
    species : str
        Species name in ``configs/species.yaml`` (must have an
        ``ensembl_dataset`` field).
    cache_dir : Path, optional
        Directory for parquet cache. Defaults to
        ``results/orthologs/cache``.
    chunked : bool, default False
        Fetch the table in per-chromosome chunks instead of one query.
        Use this for large transcriptomes (e.g. mouse) whose single-shot
        query exceeds BioMart's server-side proxy timeout. The chunks
        are concatenated and de-duplicated, so the result is identical
        to a successful single-shot fetch.
    biomart_host : str, optional
        A single BioMart ``.../biomart/martservice`` URL to query instead
        of the live mirrors. Use an Ensembl archive host (e.g.
        ``https://nov2020.archive.ensembl.org/biomart/martservice``) to
        pin the mapping to an annotation release matching an older
        deposit whose transcript IDs have since drifted out of the
        current release. When omitted, an optional
        ``tx2gene_biomart_host`` field on the species config is used,
        and finally the default live mirrors.

    Returns
    -------
    pa.Table with columns (``transcript_id``, ``gene_id``,
    ``gene_symbol``). Rows with empty transcript or gene IDs are dropped.
    """
    from src.orthologs.ensembl import _species_field

    if cache_dir is None:
        cache_dir = Path("results/orthologs/cache")
    cache_path = cache_dir / f"tx2gene_{species}.parquet"
    if cache_path.exists():
        logger.info("Loading cached tx2gene for %s from %s", species, cache_path)
        return pq.read_table(cache_path)

    dataset = _species_field(species, "ensembl_dataset")

    if biomart_host is None:
        biomart_host = _optional_species_field(species, "tx2gene_biomart_host")
    mirrors = (biomart_host,) if biomart_host else _BIOMART_MIRRORS
    if biomart_host:
        logger.info("Pinning tx2gene(%s) to archive host %s", species, biomart_host)

    if chunked:
        ensembl_species = _species_field(species, "ensembl_species")
        table, flaky = _fetch_tx2gene_chunked(dataset, ensembl_species, species, mirrors=mirrors)
    else:
        xml = _build_tx2gene_xml(dataset)
        text, flaky = _fetch_biomart_tsv(xml, species, mirrors=mirrors)
        table = _parse_tx2gene_tsv(text)

    if flaky:
        logger.warning(
            "BioMart was flaky for tx2gene(%s): %s. The successful mirror "
            "still returned data, so results are usable.",
            species,
            "; ".join(flaky),
        )

    if len(table) == 0:
        msg = f"BioMart parsed 0 tx2gene rows for {species}. Refusing to cache empty result."
        raise RuntimeError(msg)

    cache_dir.mkdir(parents=True, exist_ok=True)
    pq.write_table(table, cache_path)
    logger.info("Cached tx2gene for %s (%d rows) to %s", species, len(table), cache_path)
    return table


def _strip_version(idx) -> list[str]:
    """Strip the ``.N`` version suffix from Ensembl transcript IDs."""
    out = []
    for v in idx:
        s = str(v)
        out.append(s.split(".", 1)[0] if "." in s else s)
    return out


def aggregate_to_genes(
    adata: ad.AnnData,
    tx2gene: pa.Table,
    use_symbols: bool = False,
) -> ad.AnnData:
    """Sum-aggregate a transcript-level AnnData to gene level.

    Parameters
    ----------
    adata : ad.AnnData
        AnnData whose ``var_names`` are transcript IDs (with or without
        ``.N`` version suffix).
    tx2gene : pa.Table
        Output of :func:`fetch_tx2gene`.
    use_symbols : bool, default False
        If True, the returned AnnData uses gene symbols as ``var_names``
        (transcripts whose gene has no symbol are dropped). If False,
        Ensembl gene IDs are used and the symbol is stored in
        ``var['gene_symbol']``.

    Returns
    -------
    ad.AnnData with one row per original obs and one column per gene.
    """
    df = tx2gene.to_pandas()
    df = df.drop_duplicates("transcript_id")
    tx_to_key: dict[str, str] = {}
    tx_to_symbol: dict[str, str] = {}
    for tx, gid, sym in zip(df["transcript_id"], df["gene_id"], df["gene_symbol"], strict=True):
        if not gid:
            continue
        key = sym if (use_symbols and sym) else gid
        if use_symbols and not sym:
            continue
        tx_to_key[tx] = key
        tx_to_symbol[tx] = sym

    var_keys_in_order: list[str] = []
    seen: dict[str, int] = {}
    key_symbol: dict[str, str] = {}
    col_assignments: list[int] = []  # col index in output for each input transcript

    stripped = _strip_version(adata.var_names)
    for tx in stripped:
        key = tx_to_key.get(tx)
        if key is None:
            col_assignments.append(-1)
            continue
        if key not in seen:
            seen[key] = len(var_keys_in_order)
            var_keys_in_order.append(key)
            key_symbol[key] = tx_to_symbol.get(tx, "")
        col_assignments.append(seen[key])

    n_obs = adata.n_obs
    n_genes = len(var_keys_in_order)
    if n_genes == 0:
        msg = "No transcripts in adata mapped to any gene via tx2gene."
        raise RuntimeError(msg)

    X = adata.X
    if sp.issparse(X):
        X_csc = X.tocsc()
        out = sp.lil_matrix((n_obs, n_genes), dtype=np.float64)
        for src_col, dst_col in enumerate(col_assignments):
            if dst_col < 0:
                continue
            col = X_csc[:, src_col]
            if col.nnz == 0:
                continue
            out[:, dst_col] = out[:, dst_col] + col
        out = out.tocsr()
    else:
        X_arr = np.asarray(X, dtype=np.float64)
        out = np.zeros((n_obs, n_genes), dtype=np.float64)
        for src_col, dst_col in enumerate(col_assignments):
            if dst_col < 0:
                continue
            out[:, dst_col] += X_arr[:, src_col]

    import pandas as pd

    var_df = pd.DataFrame(index=var_keys_in_order)
    if use_symbols:
        var_df["gene_id"] = ""  # unknown after collapsing
    else:
        var_df["gene_symbol"] = [key_symbol.get(k, "") for k in var_keys_in_order]

    return ad.AnnData(
        X=out,
        obs=adata.obs.copy(),
        var=var_df,
        uns=dict(adata.uns) if adata.uns else {},
    )
