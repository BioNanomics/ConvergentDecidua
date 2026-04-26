"""Download processed data from GEO (Gene Expression Omnibus).

Strategy: processed-matrix-first. Downloads supplementary files from GEO,
preferring h5ad, MTX, or count matrices over raw FASTQs.
"""

from __future__ import annotations

import gzip
import logging
import shutil
import tarfile
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

GEO_SUPP_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/{prefix}nnn/{accession}/suppl/"
GEO_SOFT_URL = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi"

# Filename patterns we want, in priority order
_PREFERRED_EXTENSIONS = [".h5ad", ".h5", ".h5seurat", ".rds", ".loom"]
_MATRIX_PATTERNS = ["matrix.mtx", "barcodes.tsv", "features.tsv", "genes.tsv"]


def fetch_geo_dataset(
    accession: str,
    output_dir: Path,
    *,
    timeout: int = 300,
) -> list[Path]:
    """Download supplementary files for a GEO accession.

    Parameters
    ----------
    accession : str
        GEO series accession (e.g. ``"GSE127918"``).
    output_dir : Path
        Directory to save downloaded files.
    timeout : int
        HTTP request timeout in seconds.

    Returns
    -------
    list[Path]
        Paths to downloaded files.
    """
    if not accession.startswith("GSE"):
        msg = f"Expected GSE accession, got: {accession}"
        raise ValueError(msg)

    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = accession[:5]

    # Try to list supplementary files via directory listing
    supp_url = GEO_SUPP_URL.format(prefix=prefix, accession=accession)
    logger.info("Fetching supplementary file listing from %s", supp_url)

    downloaded: list[Path] = []
    try:
        file_urls = _list_supp_files(supp_url, timeout=timeout)
    except Exception:
        logger.warning("Could not list supplementary files, trying SOFT query fallback")
        file_urls = _soft_fallback(accession, timeout=timeout)

    if not file_urls:
        logger.warning("No supplementary files found for %s", accession)
        return downloaded

    # Prioritize: h5ad > h5 > MTX bundles > everything else
    prioritized = _prioritize_files(file_urls)

    for url in prioritized:
        filename = url.rsplit("/", 1)[-1]
        dest = output_dir / filename
        if dest.exists():
            logger.info("Already downloaded: %s", dest)
            downloaded.append(dest)
            continue

        logger.info("Downloading %s → %s", url, dest)
        _download_file(url, dest, timeout=timeout)
        downloaded.append(dest)

        # Extract tar.gz bundles (common for MTX format)
        if filename.endswith(".tar.gz") or filename.endswith(".tgz"):
            _extract_tar(dest, output_dir)

    return downloaded


def _list_supp_files(supp_url: str, *, timeout: int) -> list[str]:
    """Parse the FTP-style directory listing for supplementary file URLs."""
    resp = requests.get(supp_url, timeout=timeout)
    resp.raise_for_status()

    urls = []
    for line in resp.text.splitlines():
        # Simple extraction from HTML directory listing
        if 'href="' in line:
            for part in line.split('href="'):
                if part.startswith("GSE") or part.startswith("filelist"):
                    fname = part.split('"')[0]
                    if fname and not fname.endswith("/"):
                        urls.append(supp_url + fname)
    return urls


def _soft_fallback(accession: str, *, timeout: int) -> list[str]:
    """Get supplementary file URLs from SOFT metadata."""
    resp = requests.get(
        GEO_SOFT_URL,
        params={"acc": accession, "targ": "self", "form": "text", "view": "brief"},
        timeout=timeout,
    )
    resp.raise_for_status()

    urls = []
    for line in resp.text.splitlines():
        if line.startswith("!Series_supplementary_file"):
            url = line.split("=", 1)[1].strip()
            if url.startswith("http"):
                urls.append(url)
    return urls


def _prioritize_files(urls: list[str]) -> list[str]:
    """Sort URLs by preference: processed formats first."""

    def score(url: str) -> int:
        lower = url.lower()
        for i, ext in enumerate(_PREFERRED_EXTENSIONS):
            if ext in lower:
                return i
        for pattern in _MATRIX_PATTERNS:
            if pattern in lower:
                return len(_PREFERRED_EXTENSIONS)
        return len(_PREFERRED_EXTENSIONS) + 1

    return sorted(urls, key=score)


def _download_file(url: str, dest: Path, *, timeout: int) -> None:
    """Download a file with streaming."""
    resp = requests.get(url, stream=True, timeout=timeout)
    resp.raise_for_status()
    with open(dest, "wb") as fh:
        for chunk in resp.iter_content(chunk_size=8192):
            fh.write(chunk)


def _extract_tar(tar_path: Path, output_dir: Path) -> None:
    """Extract a tar.gz archive."""
    logger.info("Extracting %s", tar_path)
    with tarfile.open(tar_path, "r:gz") as tf:
        # Security: prevent path traversal
        for member in tf.getmembers():
            if member.name.startswith("/") or ".." in member.name:
                msg = f"Unsafe path in tar: {member.name}"
                raise ValueError(msg)
        tf.extractall(output_dir, filter="data")


def decompress_gz(gz_path: Path) -> Path:
    """Decompress a .gz file in place, returning the decompressed path."""
    if not gz_path.name.endswith(".gz"):
        return gz_path
    out_path = gz_path.with_suffix("")
    with gzip.open(gz_path, "rb") as f_in, open(out_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    return out_path
