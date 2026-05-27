"""Deterministic fetcher for automated_download sources.

Currently a stub: no cycle source in sources.yaml is access_class =
automated_download. This script intentionally reports skip reasons so
the workspace stays honest about provenance rather than silently
fetching nothing.
"""

from __future__ import annotations

from pathlib import Path

import yaml

HERE = Path(__file__).resolve().parent.parent
SOURCES = HERE / "sources.yaml"


def main() -> int:
    with SOURCES.open() as fh:
        manifest = yaml.safe_load(fh)

    automated = []
    skipped = []
    for section in ("cycle_sources", "contraception_sources"):
        for entry in manifest.get(section, []) or []:
            if entry.get("access_class") == "automated_download":
                automated.append(entry["id"])
            else:
                skipped.append((entry["id"], entry.get("access_class", "unknown")))

    print(f"Automated sources to fetch: {len(automated)}")
    for source_id in automated:
        print(f"  TODO fetch: {source_id}")

    print(f"\nSkipped sources ({len(skipped)}):")
    for source_id, reason in skipped:
        print(f"  {source_id}: {reason}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
