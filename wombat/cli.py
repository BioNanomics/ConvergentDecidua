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
    from pathlib import Path

    from wombat.config import load_config

    if not dataset and not fetch_all:
        console.print("[red]Specify --dataset <accession> or --all-datasets.[/red]")
        raise SystemExit(1)

    datasets = load_config("datasets")
    project_root = Path(__file__).resolve().parent.parent

    if fetch_all:
        targets = datasets
    else:
        targets = [d for d in datasets if d["accession"] == dataset]
        if not targets:
            console.print(f"[red]Unknown dataset: {dataset}[/red]")
            raise SystemExit(1)

    for ds in targets:
        acc = ds["accession"]
        raw_dir = project_root / "results" / "raw" / acc
        h5ad_path = project_root / "results" / "processed" / f"{acc}.h5ad"

        if h5ad_path.exists():
            console.print(f"[dim]Skipping {acc} — already processed[/dim]")
            continue

        console.print(f"[blue]Fetching {acc} ({ds['source']})...[/blue]")
        _fetch_one(ds, raw_dir)

        console.print(f"[blue]Converting {acc} → h5ad...[/blue]")
        from src.ingest.anndata_writer import to_anndata

        to_anndata(raw_dir, ds, h5ad_path)
        console.print(f"[green]  ✓ {h5ad_path}[/green]")


def _fetch_one(ds: dict, raw_dir: object) -> list:
    """Route download to the correct backend based on accession prefix."""
    acc = ds["accession"]
    if acc.startswith("GSE"):
        from src.ingest.geo import fetch_geo_dataset

        return fetch_geo_dataset(acc, raw_dir)
    elif acc.startswith("E-MTAB"):
        from src.ingest.arrayexpress import fetch_arrayexpress_dataset

        return fetch_arrayexpress_dataset(acc, raw_dir)
    else:
        msg = f"No downloader for accession prefix: {acc}"
        raise ValueError(msg)


@cli.command()
@click.option("--species", required=True, help="Species to QC (human or mouse).")
def qc(species: str) -> None:
    """Run quality control on processed datasets."""
    from pathlib import Path

    from wombat.config import load_config

    datasets = load_config("datasets")
    project_root = Path(__file__).resolve().parent.parent
    targets = [d for d in datasets if d["species"] == species]

    if not targets:
        console.print(f"[red]No datasets found for species: {species}[/red]")
        raise SystemExit(1)

    for ds in targets:
        acc = ds["accession"]
        h5ad_path = project_root / "results" / "processed" / f"{acc}.h5ad"
        qc_path = project_root / "results" / "qc" / f"{acc}.h5ad"

        if qc_path.exists():
            console.print(f"[dim]Skipping QC for {acc} — already done[/dim]")
            continue

        if not h5ad_path.exists():
            console.print(f"[yellow]Skipping {acc} — not yet fetched[/yellow]")
            continue

        console.print(f"[blue]Running QC on {acc} ({ds['assay']})...[/blue]")
        import anndata as ad

        adata = ad.read_h5ad(h5ad_path)
        assay = ds["assay"].lower()

        if "scrna" in assay or "snrna" in assay:
            from src.qc.scrna import qc_scrna

            adata = qc_scrna(adata, species=species)
        elif "scatac" in assay:
            from src.qc.scatac import qc_scatac

            adata = qc_scatac(adata)
        elif "bulk" in assay:
            from src.qc.bulk import qc_bulk

            adata = qc_bulk(adata)
        else:
            console.print(f"[yellow]  No QC pipeline for assay: {ds['assay']}[/yellow]")
            continue

        qc_path.parent.mkdir(parents=True, exist_ok=True)
        adata.write_h5ad(qc_path)
        console.print(f"[green]  ✓ {qc_path} ({adata.n_obs} cells)[/green]")


# ---------------------------------------------------------------------------
# Ortholog commands
# ---------------------------------------------------------------------------


@cli.group()
def orthologs() -> None:
    """Ortholog mapping commands."""


@orthologs.command("build")
@click.option("--no-gprofiler", is_flag=True, help="Skip g:Profiler cross-validation.")
def orthologs_build(no_gprofiler: bool) -> None:
    """Build ortholog backbone and orthogroup tables."""
    from pathlib import Path

    from src.orthologs.backbone import build_backbone

    project_root = Path(__file__).resolve().parent.parent
    cache_dir = project_root / "results" / "orthologs" / "cache"
    output_path = project_root / "results" / "orthologs" / "backbone.parquet"

    console.print("[blue]Building ortholog backbone (human → mouse)...[/blue]")
    backbone = build_backbone(
        source="human",
        target="mouse",
        cache_dir=cache_dir,
        output_path=output_path,
        use_gprofiler=not no_gprofiler,
    )
    console.print(f"[green]✓ Backbone: {len(backbone)} rows → {output_path}[/green]")


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
