"""Ingest path for GEO per-sample bulk-RNA dumps spanning multiple species.

Used by datasets like GSE155170 where the GEO supplement is a flat tar of
per-sample count files whose filenames carry both the species and the
biological replicate, e.g.::

    GSM4696515_Baboon_Endometrium_Individual_1.txt.gz
    GSM4696518_Hamster_Endometrium_Individual_1_Replicate_1.txt.gz
    GSM4696524_Bat_Endometrium_Individual_1.txt.gz

Because gene namespaces differ across species (each was aligned to its
own genome), we cannot concatenate samples into a single matrix without
first ortholog-mapping. This module therefore writes one ``.h5ad`` per
species (``results/processed/<accession>__<species_key>.h5ad``) and a
small "manifest" h5ad at the main ``output_path`` whose ``.uns`` records
the per-species file paths for downstream discovery.

The species mapping comes from ``dataset_meta['ingest']['filename_species_map']``,
which maps the human-readable species token in the filename to the
species key used in ``configs/species.yaml`` (e.g. ``"Bat" → "bat_carollia"``).
"""

from __future__ import annotations

import gzip
import logging
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse

logger = logging.getLogger(__name__)


# Filename: GSM<id>_<Species>_<Tissue>_Individual_<n>[_Replicate_<m>].txt(.gz)
_FILENAME_RE = re.compile(
    r"^(?P<gsm>GSM\d+)_(?P<species>[A-Za-z]+)_(?P<tissue>[A-Za-z]+)_"
    r"Individual_(?P<individual>\d+)(?:_Replicate_(?P<replicate>\d+))?"
    r"\.txt(?:\.gz)?$"
)


def write_per_sample_bulk(
    raw_dir: Path,
    dataset_meta: dict,
    output_path: Path,
) -> ad.AnnData:
    """Convert a per-sample bulk-RNA GEO dump to one h5ad per species.

    Returns
    -------
    ad.AnnData
        A tiny manifest AnnData (one obs row per per-species output file)
        written to ``output_path``. Its ``.uns['per_species_h5ads']`` maps
        species_key → file path (relative to repo root) for downstream
        QC and integration code.
    """
    accession = dataset_meta["accession"]
    ingest_cfg = dataset_meta.get("ingest") or {}
    species_map: dict[str, str] = ingest_cfg.get("filename_species_map") or {}
    if not species_map:
        msg = (
            f"{accession}: geo_per_sample_bulk requires "
            "ingest.filename_species_map in datasets.yaml"
        )
        raise ValueError(msg)

    sample_files = sorted(raw_dir.glob("GSM*.txt*"))
    if not sample_files:
        msg = f"{accession}: no GSM*.txt* sample files found in {raw_dir}"
        raise FileNotFoundError(msg)

    # Group samples by species_key
    per_species: dict[str, list[dict]] = {}
    skipped: list[str] = []
    for f in sample_files:
        m = _FILENAME_RE.match(f.name)
        if not m:
            skipped.append(f.name)
            continue
        species_token = m.group("species")
        # Normalize case for the lookup (datasets.yaml is title-case)
        species_key = species_map.get(species_token) or species_map.get(species_token.title())
        if species_key is None:
            skipped.append(f.name)
            continue
        per_species.setdefault(species_key, []).append(
            {
                "path": f,
                "sample_id": m.group("gsm"),
                "individual": m.group("individual"),
                "replicate": m.group("replicate") or "1",
                "tissue": m.group("tissue").lower(),
            }
        )

    if skipped:
        logger.warning(
            "%s: skipped %d unrecognized files (e.g. %s)",
            accession,
            len(skipped),
            skipped[:3],
        )

    per_species_h5ads: dict[str, str] = {}
    project_root = output_path.parent.parent.parent  # results/processed/.. → repo root

    for species_key, samples in per_species.items():
        logger.info(
            "%s [%s]: building h5ad from %d samples",
            accession,
            species_key,
            len(samples),
        )
        adata = _samples_to_adata(samples, accession=accession, species_key=species_key)
        adata.uns["dataset"] = {
            "accession": accession,
            "species": species_key,
            "assay": dataset_meta["assay"],
            "tissue": dataset_meta.get("tissue", ""),
            "condition": dataset_meta.get("condition", ""),
        }
        species_path = output_path.parent / f"{accession}__{species_key}.h5ad"
        species_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(species_path)
        try:
            rel = species_path.relative_to(project_root).as_posix()
        except ValueError:
            rel = str(species_path)
        per_species_h5ads[species_key] = rel
        logger.info(
            "Wrote %d samples × %d genes → %s",
            adata.n_obs,
            adata.n_vars,
            species_path,
        )

    # Manifest at output_path: one obs row per per-species file, zero genes.
    obs = pd.DataFrame(
        {
            "species": list(per_species_h5ads),
            "h5ad_path": list(per_species_h5ads.values()),
            "n_samples": [len(per_species[k]) for k in per_species_h5ads],
        },
        index=list(per_species_h5ads),
    )
    manifest = ad.AnnData(
        X=scipy.sparse.csr_matrix((len(per_species_h5ads), 0), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=pd.Index([], name="gene_id")),
    )
    manifest.uns["dataset"] = {
        "accession": accession,
        "species": dataset_meta["species"],  # "multi"
        "assay": dataset_meta["assay"],
        "tissue": dataset_meta.get("tissue", ""),
        "condition": dataset_meta.get("condition", ""),
    }
    manifest.uns["per_species_h5ads"] = per_species_h5ads
    manifest.uns["ingest_format"] = "geo_per_sample_bulk"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_h5ad(output_path)
    logger.info(
        "%s manifest: %d species → %s",
        accession,
        len(per_species_h5ads),
        output_path,
    )
    return manifest


def write_per_species_h5ad(
    raw_dir: Path,
    dataset_meta: dict,
    output_path: Path,
) -> ad.AnnData:
    """Re-index a multi-species GEO H5AD deposit as one h5ad per species.

    Used for deposits like GSE274701 where each species ships as its own
    pre-built ``.h5ad`` (e.g. ``GSE274701_Cporcellus.h5ad``). We do not
    parse the matrix at this stage — we simply copy each file into the
    project's canonical naming (``<acc>__<species_key>.h5ad``) and write
    a manifest h5ad with ``.uns['per_species_h5ads']`` pointing at them.
    Re-naming (not re-writing) keeps the original deposit pristine and
    avoids touching 100s of MB of single-cell data unnecessarily.

    ``dataset_meta['ingest']['filename_species_map']`` maps the species
    token used in the deposit filename (the part between accession and
    ``.h5ad``) to the species key from ``configs/species.yaml``.
    """
    import shutil

    accession = dataset_meta["accession"]
    ingest_cfg = dataset_meta.get("ingest") or {}
    species_map: dict[str, str] = ingest_cfg.get("filename_species_map") or {}
    if not species_map:
        msg = (
            f"{accession}: geo_per_species_h5ad requires "
            "ingest.filename_species_map in datasets.yaml"
        )
        raise ValueError(msg)

    deposit_files = sorted(raw_dir.glob(f"{accession}_*.h5ad"))
    if not deposit_files:
        msg = f"{accession}: no '{accession}_*.h5ad' deposit files found in {raw_dir}"
        raise FileNotFoundError(msg)

    per_species_h5ads: dict[str, str] = {}
    skipped: list[str] = []
    project_root = output_path.parent.parent.parent

    for f in deposit_files:
        # Token = part between "<accession>_" and ".h5ad"
        token = f.stem[len(accession) + 1 :]
        species_key = species_map.get(token) or species_map.get(token.title())
        if species_key is None:
            skipped.append(f.name)
            continue
        species_path = output_path.parent / f"{accession}__{species_key}.h5ad"
        species_path.parent.mkdir(parents=True, exist_ok=True)
        if not species_path.exists() or species_path.stat().st_size != f.stat().st_size:
            shutil.copy2(f, species_path)
        try:
            rel = species_path.relative_to(project_root).as_posix()
        except ValueError:
            rel = str(species_path)
        per_species_h5ads[species_key] = rel
        logger.info("%s [%s]: %s → %s", accession, species_key, f.name, species_path.name)

    if skipped:
        logger.warning(
            "%s: skipped %d h5ad files with unmapped species token (e.g. %s)",
            accession,
            len(skipped),
            skipped[:3],
        )

    obs = pd.DataFrame(
        {
            "species": list(per_species_h5ads),
            "h5ad_path": list(per_species_h5ads.values()),
        },
        index=list(per_species_h5ads),
    )
    manifest = ad.AnnData(
        X=scipy.sparse.csr_matrix((len(per_species_h5ads), 0), dtype=np.float32),
        obs=obs,
        var=pd.DataFrame(index=pd.Index([], name="gene_id")),
    )
    manifest.uns["dataset"] = {
        "accession": accession,
        "species": dataset_meta["species"],
        "assay": dataset_meta["assay"],
        "tissue": dataset_meta.get("tissue", ""),
        "condition": dataset_meta.get("condition", ""),
    }
    manifest.uns["per_species_h5ads"] = per_species_h5ads
    manifest.uns["ingest_format"] = "geo_per_species_h5ad"
    output_path.parent.mkdir(parents=True, exist_ok=True)
    manifest.write_h5ad(output_path)
    logger.info(
        "%s manifest: %d species → %s",
        accession,
        len(per_species_h5ads),
        output_path,
    )
    return manifest


def _samples_to_adata(
    samples: list[dict],
    *,
    accession: str,
    species_key: str,
) -> ad.AnnData:
    """Read per-sample 2-column count files and stack into samples × genes."""
    series: list[pd.Series] = []
    obs_rows: list[dict] = []
    for s in samples:
        path: Path = s["path"]
        counts = _read_two_column_counts(path)
        counts.name = s["sample_id"]
        series.append(counts)
        obs_rows.append(
            {
                "sample": s["sample_id"],
                "individual": s["individual"],
                "replicate": s["replicate"],
                "tissue": s["tissue"],
                "species": species_key,
                "dataset": accession,
            }
        )

    # Outer-join on gene IDs across samples (mostly identical for same-species
    # files, but be defensive).
    mat = pd.concat(series, axis=1).fillna(0).astype(np.float32)
    # samples × genes
    X = scipy.sparse.csr_matrix(mat.T.values)
    obs = pd.DataFrame(obs_rows, index=[s["sample_id"] for s in samples])
    var = pd.DataFrame(index=mat.index)
    var.index.name = "gene_id"
    return ad.AnnData(X=X, obs=obs, var=var)


def _read_two_column_counts(path: Path) -> pd.Series:
    """Read a per-sample count file as a Series indexed by gene/transcript id.

    Supports two on-disk schemas:
      1. Plain ``gene_id<TAB>count`` (with or without header).
      2. kallisto ``abundance.tsv``-style: 6 columns
         ``target_id, target_id version, length, eff_length, est_counts, tpm``.
         The ``target_id`` is used as the index and ``est_counts`` as the
         value. (Used by GSE155170; ids are transcript-level — gene rollup
         is the QC layer's job.)
    """
    opener = gzip.open if path.name.endswith(".gz") else open
    with opener(path, "rt") as fh:
        df = pd.read_csv(fh, sep="\t")

    if "target_id" in df.columns and "est_counts" in df.columns:
        idx = df["target_id"].astype(str)
        vals = pd.to_numeric(df["est_counts"], errors="coerce").fillna(0.0)
    elif df.shape[1] >= 2:
        # Headerless or plain 2-col: first = id, second = count
        idx = df.iloc[:, 0].astype(str)
        vals = pd.to_numeric(df.iloc[:, 1], errors="coerce").fillna(0.0)
    else:
        msg = f"Cannot parse count file {path}: shape={df.shape}"
        raise ValueError(msg)

    s = pd.Series(vals.values, index=idx, dtype=np.float32)
    s = s[~s.index.duplicated(keep="first")]
    return s
