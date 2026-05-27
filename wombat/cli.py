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


@orthologs.command("synteny-check")
@click.option(
    "--symbols",
    default="",
    help="Comma-separated anchor-species gene symbols. Defaults to "
    "configs/markers.yaml::protected_core.",
)
@click.option(
    "--targets",
    default="mouse",
    help="Comma-separated target species short names (e.g. mouse,rat).",
)
@click.option("--anchor", default="human", help="Anchor species short name.")
def orthologs_synteny_check(symbols: str, targets: str, anchor: str) -> None:
    """Per-locus 1:1 ortholog + alignment % identity via Ensembl REST.

    Replaces the manual CGV eyeball loop for the protected-core
    decidual panel. See ``src/orthologs/synteny.py``.
    """
    from pathlib import Path

    from src.orthologs.synteny import run_synteny_check

    project_root = Path(__file__).resolve().parent.parent
    output_path = project_root / "results" / "orthologs" / "synteny_at_core_loci.parquet"

    sym_list = [s.strip() for s in symbols.split(",") if s.strip()] or None
    target_list = [t.strip() for t in targets.split(",") if t.strip()]

    console.print(
        f"[blue]Synteny check: anchor={anchor}, targets={target_list}, "
        f"symbols={'protected_core' if sym_list is None else sym_list}[/blue]"
    )
    table = run_synteny_check(
        output_path,
        symbols=sym_list,
        target_species=target_list,
        anchor_species=anchor,
    )
    df = table.to_pandas()
    hits = int(df["alignment_present"].sum())
    console.print(f"[green]✓ {hits}/{len(df)} alignment-present rows → {output_path}[/green]")


# ---------------------------------------------------------------------------
# Integration & scoring commands
# ---------------------------------------------------------------------------


@cli.command()
@click.option("--mode", default="stromal", help="Integration mode (e.g. stromal).")
@click.option("--method", default="harmony", help="Integration method: harmony or scvi.")
@click.option(
    "--orthology-tier",
    type=click.Choice(["1", "12"]),
    default="1",
    help="Ortholog tier: 1 = high-confidence 1:1 only; 12 = include Tier 2 orthogroups.",
)
def integrate(mode: str, method: str, orthology_tier: str) -> None:
    """Integrate datasets across species."""
    from pathlib import Path

    import anndata as ad

    from src.cell_states.annotate import annotate_cell_types
    from src.cell_states.integrate import integrate_stromal
    from src.cell_states.subset import subset_stromal
    from wombat.config import load_config

    project_root = Path(__file__).resolve().parent.parent
    backbone_path = project_root / "results" / "orthologs" / "backbone.parquet"
    qc_dir = project_root / "results" / "qc"
    canonical_path = project_root / "results" / "integrated" / f"{mode}_cross_species.h5ad"
    legacy_alias = project_root / "results" / "integrated" / f"{mode}_{method}.h5ad"

    if not backbone_path.exists():
        console.print("[red]Backbone not found — run 'wombat orthologs build' first[/red]")
        raise SystemExit(1)

    datasets = load_config("datasets")
    stromal_list = []

    for ds in datasets:
        acc = ds["accession"]
        qc_path = qc_dir / f"{acc}.h5ad"
        if not qc_path.exists():
            console.print(f"[yellow]Skipping {acc} — not QC'd[/yellow]")
            continue

        # Only scRNA-seq datasets for cell-state integration
        if "scrna" not in ds["assay"].lower() and "snrna" not in ds["assay"].lower():
            console.print(f"[dim]Skipping {acc} — not scRNA-seq[/dim]")
            continue

        console.print(f"[blue]Annotating {acc}...[/blue]")
        adata = ad.read_h5ad(qc_path)
        adata = annotate_cell_types(
            adata,
            species=ds["species"],
            backbone_path=str(backbone_path),
        )

        if mode == "stromal":
            adata = subset_stromal(adata)

        if adata.n_obs > 0:
            stromal_list.append(adata)
            console.print(f"[green]  ✓ {acc}: {adata.n_obs} {mode} cells[/green]")

    if not stromal_list:
        console.print("[red]No datasets with qualifying cells[/red]")
        raise SystemExit(1)

    console.print(f"[blue]Integrating {len(stromal_list)} datasets ({method})...[/blue]")
    markers_cfg = load_config("markers")
    protected_core = markers_cfg.get("protected_core") if isinstance(markers_cfg, dict) else None
    if protected_core:
        console.print(f"[dim]  Protected-core panel: {protected_core}[/dim]")
    integrated = integrate_stromal(
        stromal_list,
        backbone_path,
        method=method,
        orthology_tier=int(orthology_tier),
        protected_core=protected_core,
    )

    canonical_path.parent.mkdir(parents=True, exist_ok=True)
    integrated.write_h5ad(canonical_path)
    console.print(f"[green]✓ Integrated: {integrated.n_obs} cells → {canonical_path}[/green]")

    # Back-compat alias for the legacy {mode}_{method}.h5ad name.
    # Use a symlink when possible; fall back to a copy on filesystems
    # that don't support it (or where one already exists as a real file).
    try:
        if legacy_alias.is_symlink() or legacy_alias.exists():
            legacy_alias.unlink()
        legacy_alias.symlink_to(canonical_path.name)
        console.print(f"[dim]  ↳ alias: {legacy_alias.name} -> {canonical_path.name}[/dim]")
    except OSError:
        import shutil

        shutil.copy2(canonical_path, legacy_alias)
        console.print(f"[dim]  ↳ alias copy: {legacy_alias.name}[/dim]")


@cli.command("score-decidua")
def score_decidua() -> None:
    """Compute decidualization scores across all 8 modules."""
    from pathlib import Path

    import anndata as ad

    from src.scoring.engine import score_all_modules
    from src.scoring.gene_sets import load_score_gene_sets
    from src.scoring.reports import generate_score_report

    project_root = Path(__file__).resolve().parent.parent
    backbone_path = project_root / "results" / "orthologs" / "backbone.parquet"
    integrated_path = project_root / "results" / "integrated" / "stromal_harmony.h5ad"
    report_dir = project_root / "results" / "reports" / "scoring"

    if not integrated_path.exists():
        console.print("[red]Integrated data not found — run 'wombat integrate' first[/red]")
        raise SystemExit(1)

    console.print("[blue]Loading integrated data...[/blue]")
    adata = ad.read_h5ad(integrated_path)
    gene_sets = load_score_gene_sets()

    console.print(f"[blue]Scoring {len(gene_sets)} modules...[/blue]")
    species = adata.obs["species"].iloc[0] if "species" in adata.obs.columns else "human"
    bp = backbone_path if backbone_path.exists() else None
    adata = score_all_modules(adata, gene_sets, species=species, backbone_path=bp)

    # Save scored object
    scored_path = project_root / "results" / "scored" / "stromal_scored.h5ad"
    scored_path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(scored_path)
    console.print(f"[green]✓ Scored: {scored_path}[/green]")

    # Generate reports
    console.print("[blue]Generating scoring report...[/blue]")
    score_cols = list(gene_sets.keys())
    generate_score_report(adata, score_cols, report_dir)
    console.print(f"[green]✓ Report: {report_dir}[/green]")


# ---------------------------------------------------------------------------
# Reports command
# ---------------------------------------------------------------------------


@cli.command("generate-reports")
def generate_reports() -> None:
    """Generate all pipeline reports and release manifest."""
    from pathlib import Path

    from src.reports.coverage import generate_coverage_report
    from src.reports.integration_qc import generate_integration_qc
    from src.reports.manifest import generate_manifest
    from src.reports.methods import generate_methods_report
    from src.reports.ortholog_report import generate_ortholog_report
    from src.reports.qc_report import generate_qc_report

    project_root = Path(__file__).resolve().parent.parent
    results_dir = project_root / "results"
    report_dir = results_dir / "reports"

    console.print("[blue]Generating reports...[/blue]")

    generate_methods_report(report_dir / "methods.md")
    console.print("[green]  ✓ Methods report[/green]")

    generate_coverage_report(results_dir, report_dir / "coverage.md")
    console.print("[green]  ✓ Coverage report[/green]")

    generate_qc_report(results_dir, report_dir / "qc_summary.md")
    console.print("[green]  ✓ QC report[/green]")

    generate_ortholog_report(results_dir, report_dir / "orthologs.md")
    console.print("[green]  ✓ Ortholog report[/green]")

    generate_integration_qc(
        results_dir / "integrated" / "stromal_cross_species.h5ad",
        report_dir / "integration_qc.md",
        backbone_path=results_dir / "orthologs" / "backbone.parquet",
        processed_dir=results_dir / "processed",
    )
    console.print("[green]  ✓ Integration QC report[/green]")

    generate_manifest(results_dir, report_dir / "manifest.md")
    console.print("[green]  ✓ Release manifest[/green]")

    console.print(f"[green]All reports → {report_dir}/[/green]")


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
