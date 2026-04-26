"""Wombat CLI — command-line interface for ConvergentDecidua.

Entry point: ``wombat``
"""

from __future__ import annotations

import click
from rich.console import Console

from wombat.config import validate_all

console = Console()


@click.group()
@click.option("-v", "--verbose", count=True, help="Increase verbosity (-v, -vv).")
def cli(verbose: int) -> None:
    """Wombat — ConvergentDecidua workflow CLI."""
    from wombat.logging import setup_logging

    setup_logging(verbose)


# ---------------------------------------------------------------------------
# Config commands
# ---------------------------------------------------------------------------


@cli.command()
def init() -> None:
    """Validate that all required configs exist."""
    from pathlib import Path

    configs_dir = Path(__file__).resolve().parent.parent / "configs"
    required = ["datasets.yaml", "species.yaml", "markers.yaml"]
    missing = [f for f in required if not (configs_dir / f).exists()]
    if missing:
        console.print(f"[red]Missing configs: {missing}[/red]")
        raise SystemExit(1)
    console.print("[green]All required configs present.[/green]")


@cli.command("validate-config")
def validate_config() -> None:
    """Validate all YAML configuration files."""
    errors = validate_all()
    if errors:
        for err in errors:
            console.print(f"[red]  ✗ {err}[/red]")
        raise SystemExit(1)
    console.print("[green]All configs valid.[/green]")


@cli.command("build-registry")
def build_registry() -> None:
    """Export dataset registry to Parquet and CSV."""
    from pathlib import Path

    import pyarrow as pa
    import pyarrow.parquet as pq

    from wombat.config import load_config

    datasets = load_config("datasets")
    table = pa.Table.from_pylist(datasets)

    out_dir = Path(__file__).resolve().parent.parent / "results"
    out_dir.mkdir(parents=True, exist_ok=True)

    pq.write_table(table, out_dir / "registry.parquet")
    table.to_pandas().to_csv(out_dir / "registry.csv", index=False)
    console.print(
        f"[green]Registry exported: {len(datasets)} datasets → results/registry.parquet, results/registry.csv[/green]"
    )


# ---------------------------------------------------------------------------
# Data commands
# ---------------------------------------------------------------------------


@cli.command()
@click.option("--dataset", help="Accession to fetch (e.g. GSE127918).")
@click.option("--all-datasets", "fetch_all", is_flag=True, help="Fetch all priority datasets.")
def fetch(dataset: str | None, fetch_all: bool) -> None:
    """Download datasets and convert to standardized AnnData."""
    if not dataset and not fetch_all:
        console.print("[red]Specify --dataset <accession> or --all-datasets.[/red]")
        raise SystemExit(1)
    console.print(f"[yellow]fetch: not yet implemented (dataset={dataset})[/yellow]")


@cli.command()
@click.option("--species", required=True, help="Species to QC (human or mouse).")
def qc(species: str) -> None:
    """Run quality control on processed datasets."""
    console.print(f"[yellow]qc: not yet implemented (species={species})[/yellow]")


# ---------------------------------------------------------------------------
# Ortholog commands
# ---------------------------------------------------------------------------


@cli.group()
def orthologs() -> None:
    """Ortholog mapping commands."""


@orthologs.command("build")
def orthologs_build() -> None:
    """Build ortholog backbone and orthogroup tables."""
    console.print("[yellow]orthologs build: not yet implemented[/yellow]")


# ---------------------------------------------------------------------------
# Integration & scoring commands
# ---------------------------------------------------------------------------


@cli.command()
@click.option("--mode", default="stromal", help="Integration mode (e.g. stromal).")
def integrate(mode: str) -> None:
    """Integrate datasets across species."""
    console.print(f"[yellow]integrate: not yet implemented (mode={mode})[/yellow]")


@cli.command("score-decidua")
def score_decidua() -> None:
    """Compute decidualization scores across all modules."""
    console.print("[yellow]score-decidua: not yet implemented[/yellow]")


# ---------------------------------------------------------------------------
# Atlas command
# ---------------------------------------------------------------------------


@cli.command("serve-atlas")
@click.option("--port", default=8501, help="Streamlit port.")
def serve_atlas(port: int) -> None:
    """Launch the DecidualAtlas Streamlit app."""
    import subprocess
    import sys
    from pathlib import Path

    app_path = Path(__file__).resolve().parent.parent / "decidual_atlas" / "app.py"
    if not app_path.exists():
        console.print("[red]Atlas app not found.[/red]")
        raise SystemExit(1)
    subprocess.run(
        [sys.executable, "-m", "streamlit", "run", str(app_path), "--server.port", str(port)],
        check=True,
    )
