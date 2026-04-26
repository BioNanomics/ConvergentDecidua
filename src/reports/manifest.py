"""Release manifest with checksums.

Generates a manifest of all result files with SHA-256 checksums
for reproducibility and data integrity verification.
"""

from __future__ import annotations

import hashlib
import logging
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# File patterns to include in the manifest
RESULT_PATTERNS = [
    "**/*.h5ad",
    "**/*.parquet",
    "**/*.csv",
    "**/*.png",
    "**/*.md",
]


def generate_manifest(results_dir: Path, output_path: Path) -> pd.DataFrame:
    """Generate release manifest with checksums.

    Parameters
    ----------
    results_dir : Path
        Path to the results directory.
    output_path : Path
        Where to write the manifest.

    Returns
    -------
    pd.DataFrame
        Manifest with file paths, sizes, and checksums.
    """
    records = []

    for pattern in RESULT_PATTERNS:
        for path in sorted(results_dir.glob(pattern)):
            if path.is_file():
                stat = path.stat()
                records.append(
                    {
                        "path": str(path.relative_to(results_dir)),
                        "size_bytes": stat.st_size,
                        "size_human": _human_size(stat.st_size),
                        "sha256": _sha256(path),
                        "modified": datetime.fromtimestamp(
                            stat.st_mtime, tz=timezone.utc
                        ).isoformat(),
                    }
                )

    df = pd.DataFrame(records)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        ts = datetime.now(tz=timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
        fh.write(f"# Release Manifest\n\n*Generated: {ts}*\n\n")
        fh.write(f"**Total files**: {len(df)}\n\n")
        if len(df) > 0:
            total_bytes = df["size_bytes"].sum()
            fh.write(f"**Total size**: {_human_size(total_bytes)}\n\n")
            fh.write(df[["path", "size_human", "sha256"]].to_markdown(index=False))
        else:
            fh.write("No result files found.\n")
        fh.write("\n")

    # Also write CSV for programmatic use
    csv_path = output_path.with_suffix(".csv")
    df.to_csv(csv_path, index=False)

    logger.info("Manifest: %d files → %s", len(df), output_path)
    return df


def _sha256(path: Path) -> str:
    """Compute SHA-256 hash of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def _human_size(n: int) -> str:
    """Convert bytes to human-readable string."""
    for unit in ("B", "KB", "MB", "GB"):
        if n < 1024:
            return f"{n:.1f} {unit}"
        n /= 1024
    return f"{n:.1f} TB"
