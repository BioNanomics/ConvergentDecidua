"""Cross-validate orthologs using g:Profiler g:Orth API.

Provides an independent validation of Ensembl-derived orthologs
and fills gaps where BioMart data is incomplete.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

logger = logging.getLogger(__name__)

GORTH_URL = "https://biit.cs.ut.ee/gprofiler/api/orth/orth/"

_ORGANISM_MAP = {
    "human": "hsapiens",
    "mouse": "mmusculus",
}


def fetch_gprofiler_orthologs(
    gene_list: list[str],
    source: str = "human",
    target: str = "mouse",
    cache_dir: Path | None = None,
) -> pa.Table:
    """Query g:Profiler g:Orth for ortholog mapping.

    Parameters
    ----------
    gene_list : list[str]
        List of source gene symbols or Ensembl IDs.
    source : str
        Source species name.
    target : str
        Target species name.
    cache_dir : Path, optional
        Directory to cache results as Parquet.

    Returns
    -------
    pa.Table
        Table with source_gene, target_gene, target_symbol, mapping_type.
    """
    import requests

    src_org = _ORGANISM_MAP.get(source)
    tgt_org = _ORGANISM_MAP.get(target)
    if not src_org or not tgt_org:
        msg = f"Unknown species pair: {source} → {target}"
        raise ValueError(msg)

    # Check cache
    if cache_dir is not None:
        cache_path = cache_dir / f"gprofiler_{source}_{target}.parquet"
        if cache_path.exists():
            logger.info("Loading cached g:Profiler orthologs from %s", cache_path)
            return pq.read_table(cache_path)

    # Query in batches (API limit ~1000 genes per request)
    batch_size = 500
    all_results = []

    for i in range(0, len(gene_list), batch_size):
        batch = gene_list[i : i + batch_size]
        payload = {
            "organism": src_org,
            "target": tgt_org,
            "query": batch,
        }

        resp = requests.post(GORTH_URL, json=payload, timeout=120)
        resp.raise_for_status()
        data = resp.json()

        for item in data.get("result", []):
            for mapping in item.get("ortholog_ensg", []):
                all_results.append(
                    {
                        "source_gene": item.get("incoming", ""),
                        "target_gene": mapping,
                        "target_symbol": "",  # g:Orth returns Ensembl IDs
                        "mapping_type": "gprofiler",
                    }
                )

        logger.debug("g:Orth batch %d–%d: %d mappings", i, i + len(batch), len(all_results))

    table = pa.table(
        {
            "source_gene": [r["source_gene"] for r in all_results],
            "target_gene": [r["target_gene"] for r in all_results],
            "target_symbol": [r["target_symbol"] for r in all_results],
            "mapping_type": [r["mapping_type"] for r in all_results],
        }
    )
    logger.info("g:Profiler: %d ortholog mappings", len(table))

    # Cache
    if cache_dir is not None:
        cache_dir.mkdir(parents=True, exist_ok=True)
        pq.write_table(table, cache_path)
        logger.info("Cached to %s", cache_path)

    return table
