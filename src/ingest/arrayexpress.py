"""Download processed data from ArrayExpress / BioStudies.

Handles E-MTAB accessions via the BioStudies REST API.
"""

from __future__ import annotations

import logging
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

BIOSTUDIES_API = "https://www.ebi.ac.uk/biostudies/api/v1/studies/{accession}"
BIOSTUDIES_FILES = "https://www.ebi.ac.uk/biostudies/files/{accession}/{path}"


def fetch_arrayexpress_dataset(
    accession: str,
    output_dir: Path,
    *,
    timeout: int = 300,
) -> list[Path]:
    """Download supplementary/processed files for an ArrayExpress accession.

    Parameters
    ----------
    accession : str
        ArrayExpress accession (e.g. ``"E-MTAB-11491"``).
    output_dir : Path
        Directory to save downloaded files.
    timeout : int
        HTTP request timeout in seconds.

    Returns
    -------
    list[Path]
        Paths to downloaded files.
    """
    if not accession.startswith("E-MTAB"):
        msg = f"Expected E-MTAB accession, got: {accession}"
        raise ValueError(msg)

    output_dir.mkdir(parents=True, exist_ok=True)

    # Get study metadata from BioStudies API
    logger.info("Querying BioStudies for %s", accession)
    resp = requests.get(
        BIOSTUDIES_API.format(accession=accession),
        timeout=timeout,
    )
    resp.raise_for_status()
    study = resp.json()

    # Extract file paths from study sections
    file_paths = _extract_file_paths(study)
    logger.info("Found %d files for %s", len(file_paths), accession)

    downloaded: list[Path] = []
    for fpath in file_paths:
        filename = fpath.rsplit("/", 1)[-1]
        dest = output_dir / filename
        if dest.exists():
            logger.info("Already downloaded: %s", dest)
            downloaded.append(dest)
            continue

        url = BIOSTUDIES_FILES.format(accession=accession, path=fpath)
        logger.info("Downloading %s → %s", url, dest)
        _download_file(url, dest, timeout=timeout)
        downloaded.append(dest)

    return downloaded


def _extract_file_paths(study: dict) -> list[str]:
    """Extract downloadable file paths from a BioStudies JSON response."""
    paths: list[str] = []

    def _walk(section: dict) -> None:
        # Check for file references in this section
        if "files" in section:
            for file_entry in section["files"]:
                if isinstance(file_entry, list):
                    for f in file_entry:
                        if isinstance(f, dict) and "path" in f:
                            paths.append(f["path"])
                elif isinstance(file_entry, dict) and "path" in file_entry:
                    paths.append(file_entry["path"])
        # Recurse into subsections
        if "subsections" in section:
            for sub in section["subsections"]:
                if isinstance(sub, list):
                    for s in sub:
                        if isinstance(s, dict):
                            _walk(s)
                elif isinstance(sub, dict):
                    _walk(sub)

    if "section" in study:
        _walk(study["section"])

    return paths


def _download_file(url: str, dest: Path, *, timeout: int) -> None:
    """Download a file with streaming."""
    resp = requests.get(url, stream=True, timeout=timeout)
    resp.raise_for_status()
    with open(dest, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=8192):
            fh.write(chunk)
