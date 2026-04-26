"""Ortholog mapping report.

Summarizes the ortholog backbone: tier counts, coverage, key markers.
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)

# Key markers that should be in the backbone (README §8-9)
KEY_MARKERS = [
    "PRL",
    "IGFBP1",
    "PGR",
    "FOXO1",
    "HOXA10",
    "HAND2",
    "WNT4",
    "BMP2",
    "ESR1",
    "VEGFA",
    "MMP2",
    "IL15",
]


def generate_ortholog_report(results_dir: Path, output_path: Path) -> pd.DataFrame | None:
    """Generate ortholog mapping report.

    Parameters
    ----------
    results_dir : Path
        Path to the results directory.
    output_path : Path
        Where to write the report.

    Returns
    -------
    pd.DataFrame or None
        Backbone summary, or None if backbone doesn't exist.
    """
    backbone_path = results_dir / "orthologs" / "backbone.parquet"
    if not backbone_path.exists():
        logger.warning("Backbone not found at %s", backbone_path)
        return None

    df = pd.read_parquet(backbone_path)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as fh:
        fh.write("# Ortholog Mapping Report\n\n")

        # Summary stats
        fh.write("## Summary\n\n")
        fh.write(f"- **Total mappings**: {len(df)}\n")
        for tier in sorted(df["tier"].unique()):
            n = (df["tier"] == tier).sum()
            fh.write(f"- **Tier {tier}**: {n} mappings\n")
        fh.write(f"- **Unique source genes**: {df['source_symbol'].nunique()}\n")
        fh.write(f"- **Unique target genes**: {df['target_symbol'].nunique()}\n")

        if "gprofiler_confirmed" in df.columns:
            confirmed = df[df["tier"] == 1]["gprofiler_confirmed"].sum()
            total_t1 = (df["tier"] == 1).sum()
            fh.write(f"- **g:Profiler confirmed (Tier 1)**: {confirmed}/{total_t1}\n")

        # Key marker check
        fh.write("\n## Key Marker Coverage\n\n")
        fh.write("| Marker | In Backbone | Tier | Target Symbol |\n")
        fh.write("|--------|-------------|------|---------------|\n")
        for marker in KEY_MARKERS:
            row = df[df["source_symbol"] == marker]
            if len(row) > 0:
                best = row.sort_values("tier").iloc[0]
                fh.write(f"| {marker} | ✅ | {best['tier']} | {best['target_symbol']} |\n")
            else:
                fh.write(f"| {marker} | ❌ | — | — |\n")

        # Orthology type distribution
        fh.write("\n## Orthology Types\n\n")
        type_counts = df["orthology_type"].value_counts()
        fh.write(type_counts.to_markdown())
        fh.write("\n")

    logger.info("Ortholog report → %s", output_path)
    return df
