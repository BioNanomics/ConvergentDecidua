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


@cli.command("score-bulk")
@click.option(
    "--accession",
    help="Single bulk accession to score (default: all bulk datasets in configs/datasets.yaml).",
)
def score_bulk(accession: str | None) -> None:
    """Score bulk RNA-seq datasets and report monotonicity vs. time (Q3.1)."""
    from pathlib import Path

    import anndata as ad

    from src.scoring.bulk import monotonicity
    from src.scoring.bulk import score_bulk as _score_bulk
    from src.scoring.gene_sets import load_score_gene_sets
    from src.scoring.reports import generate_bulk_score_report
    from wombat.config import load_config

    project_root = Path(__file__).resolve().parent.parent
    backbone_path = project_root / "results" / "orthologs" / "backbone.parquet"
    report_dir = project_root / "results" / "reports" / "scoring"
    scored_dir = project_root / "results" / "scored"
    scored_dir.mkdir(parents=True, exist_ok=True)

    datasets = load_config("datasets")
    bulk = [d for d in datasets if "bulk" in d["assay"].lower()]
    if accession:
        bulk = [d for d in bulk if d["accession"] == accession]
    if not bulk:
        console.print("[red]No bulk datasets matched.[/red]")
        raise SystemExit(1)

    gene_sets = load_score_gene_sets()
    bp = backbone_path if backbone_path.exists() else None

    scored: dict[str, ad.AnnData] = {}
    tables = {}
    for ds in bulk:
        acc = ds["accession"]
        qc_path = project_root / "results" / "qc" / f"{acc}.h5ad"
        if not qc_path.exists():
            console.print(f"[yellow]Skip {acc} — no QC h5ad at {qc_path}[/yellow]")
            continue
        console.print(f"[blue]Scoring {acc} ({ds['species']})...[/blue]")
        adata = ad.read_h5ad(qc_path)
        adata = _score_bulk(adata, species=ds["species"], backbone_path=bp, gene_sets=gene_sets)
        out = scored_dir / f"{acc}_bulk_scored.h5ad"
        adata.write_h5ad(out)
        console.print(f"[green]  ✓ {out}[/green]")

        table = monotonicity(adata, list(gene_sets.keys()))
        tables[acc] = table
        scored[acc] = adata
        dec = table.loc["decidual_score"] if "decidual_score" in table.index else None
        if dec is not None:
            console.print(
                f"[dim]    decidual_score: rho={dec['rho']:+.3f}, "
                f"pval={dec['pval']:.3g}, monotonic={dec['monotonic']}[/dim]"
            )

    if scored:
        generate_bulk_score_report(scored, tables, report_dir)
        console.print(
            f"[green]✓ Bulk scoring report → {report_dir}/bulk_scoring_report.md[/green]"
        )


@cli.command("score-null")
@click.option("--n-permutations", default=200, show_default=True, type=int)
@click.option("--seed", default=0, show_default=True, type=int)
@click.option(
    "--group-key",
    default="cell_type",
    show_default=True,
    help="``.obs`` column to aggregate within (e.g. cell_type, cell_state).",
)
def score_null(n_permutations: int, seed: int, group_key: str) -> None:
    """Permutation null + per-(species, cell_state) FDR (Q3.2)."""
    from pathlib import Path

    import anndata as ad

    from src.scoring.gene_sets import load_score_gene_sets
    from src.scoring.null import NullConfig, score_with_null

    project_root = Path(__file__).resolve().parent.parent
    scored_path = project_root / "results" / "scored" / "stromal_scored.h5ad"
    if not scored_path.exists():
        console.print(f"[red]Missing {scored_path}. Run `wombat score-decidua` first.[/red]")
        raise SystemExit(1)

    backbone_path = project_root / "results" / "orthologs" / "backbone.parquet"
    # The integrated atlas is symbol-harmonised to human (mouse cells share
    # the same human symbol var-namespace), so no per-species mapping is
    # needed here. Bulk/per-dataset scoring uses the backbone elsewhere.
    _ = backbone_path  # noqa: F841 — kept for future per-dataset (non-integrated) callers

    console.print(f"[blue]Loading {scored_path}...[/blue]")
    adata = ad.read_h5ad(scored_path)
    gene_sets = load_score_gene_sets()
    species_to_backbone: dict[str, str | None] = dict.fromkeys(
        adata.obs["species"].astype(str).unique(), None
    )

    console.print(
        f"[blue]Permutation null: {n_permutations} draws × "
        f"{len(gene_sets)} modules × {adata.obs['species'].nunique()} species[/blue]"
    )
    table = score_with_null(
        adata,
        gene_sets,
        species_to_backbone=species_to_backbone,
        config=NullConfig(n_permutations=n_permutations, seed=seed, group_key=group_key),
    )

    out_dir = project_root / "results" / "reports" / "scoring"
    out_dir.mkdir(parents=True, exist_ok=True)
    csv_path = out_dir / "permutation_fdr.csv"
    md_path = out_dir / "permutation_fdr.md"
    table.to_csv(csv_path, index=False)

    with open(md_path, "w") as fh:
        fh.write("# Permutation-null FDR (Q3.2)\n\n")
        fh.write(
            f"- n_permutations = **{n_permutations}**, seed = {seed}\n"
            "- size-matched random gene sets drawn from "
            "``adata.var_names`` per species\n"
            "- one-sided absolute-deviation test, BH-corrected across all rows\n\n"
        )
        sig = table[table["fdr"] < 0.05]
        fh.write(f"## Summary\n\n- rows tested: {len(table)}\n")
        fh.write(f"- rows with FDR < 0.05: **{len(sig)}**\n")
        conserved = (
            sig.groupby("module")["species"]
            .nunique()
            .pipe(lambda s: s[s == adata.obs["species"].nunique()])
            .index.tolist()
        )
        fh.write(
            f"- modules significant in **all** species (Q3.3 conserved-pool seed): "
            f"{conserved or '∅'}\n\n"
        )
        fh.write("## Full table\n\n")
        fh.write(table.to_markdown(index=False))
    console.print(f"[green]✓ {csv_path}[/green]")
    console.print(f"[green]✓ {md_path}[/green]")


@cli.command("classify-conservation")
@click.option("--fdr-threshold", default=0.05, show_default=True, type=float)
def classify_conservation_cmd(fdr_threshold: float) -> None:
    """Classify modules as conserved / biased / divergent (Q3.3)."""
    from pathlib import Path

    import pandas as pd

    from src.scoring.conservation import classify_conservation, summarise_modules

    project_root = Path(__file__).resolve().parent.parent
    fdr_csv = project_root / "results" / "reports" / "scoring" / "permutation_fdr.csv"
    if not fdr_csv.exists():
        console.print(f"[red]Missing {fdr_csv}. Run `wombat score-null` first.[/red]")
        raise SystemExit(1)

    fdr_table = pd.read_csv(fdr_csv)
    detail = classify_conservation(fdr_table, fdr_threshold=fdr_threshold)
    summary = summarise_modules(detail)

    out_dir = project_root / "results" / "reports"
    detail_csv = out_dir / "conservation_table.csv"
    summary_csv = out_dir / "conservation_summary.csv"
    md_path = out_dir / "conservation_table.md"

    detail.to_csv(detail_csv, index=False)
    summary.to_csv(summary_csv, index=False)

    with open(md_path, "w") as fh:
        fh.write("# Conservation table (Q3.3)\n\n")
        fh.write(
            f"Classes derived from `permutation_fdr.csv` at "
            f"FDR < **{fdr_threshold}** per (module × cell group × species).\n\n"
            "Classes: `conserved-{up,down}` (sig in all species, same sign),\n"
            "`divergent` (sig in all species, opposite signs),\n"
            "`{species}-biased-{up,down}`, or `neutral`.\n\n"
        )
        fh.write("## Module-level summary\n\n")
        fh.write(summary.to_markdown(index=False))
        fh.write("\n\n## Per-group detail\n\n")
        fh.write(detail.to_markdown(index=False))

    console.print(f"[green]✓ {detail_csv}[/green]")
    console.print(f"[green]✓ {summary_csv}[/green]")
    console.print(f"[green]✓ {md_path}[/green]")
    n_conserved = int(summary["summary_class"].str.startswith("conserved-").sum())
    console.print(f"[bold green]Conserved modules: {n_conserved}/{len(summary)}[/bold green]")


@cli.command("rank-regulators")
@click.option("--cap", default=25, show_default=True, type=int)
@click.option("--tf-list", default=None, type=str, help="Path to TF symbol list.")
def rank_regulators_cmd(cap: int, tf_list: str | None) -> None:
    """Rank candidate decidualization regulators (Q3.4)."""
    from pathlib import Path

    import anndata as ad

    from src.cell_states.regulators import (
        RegulatorConfig,
        load_tf_list,
        rank_regulators,
        split_regulators,
    )

    project_root = Path(__file__).resolve().parent.parent
    scored_path = project_root / "results" / "scored" / "stromal_scored.h5ad"
    if not scored_path.exists():
        console.print(f"[red]Missing {scored_path}. Run `wombat score-decidua` first.[/red]")
        raise SystemExit(1)

    console.print(f"[blue]Loading {scored_path}...[/blue]")
    adata = ad.read_h5ad(scored_path)
    tfs = load_tf_list(tf_list)
    console.print(f"[blue]Ranking {len(tfs)} candidate TFs vs decidual_score...[/blue]")
    ranked = rank_regulators(adata, tfs, config=RegulatorConfig(cap=cap))
    splits = split_regulators(ranked, cap=cap)

    out_dir = project_root / "results" / "reports"
    full_csv = out_dir / "regulators_full_ranking.csv"
    ranked.to_csv(full_csv, index=False)
    for name, df in splits.items():
        df.to_csv(out_dir / f"regulators_{name}.csv", index=False)

    md_path = out_dir / "regulators.md"
    with open(md_path, "w") as fh:
        fh.write("# Candidate decidualization regulators (Q3.4)\n\n")
        fh.write(
            "Spearman correlation between each Lambert-2018 human TF's\n"
            "expression and ``decidual_score`` within the decidual lineage\n"
            "(``pre_decidual_stromal`` + ``decidual_stromal`` + "
            "``senescent_decidual``), per species. Ranking is by |rho|.\n\n"
        )
        for name, df in splits.items():
            fh.write(f"## {name.replace('_', ' ').title()} (top {len(df)})\n\n")
            cols = [
                c
                for c in [
                    "tf",
                    "human_rho",
                    "mouse_rho",
                    "human_rank",
                    "mouse_rank",
                    "mean_rank",
                    "rank_gap",
                ]
                if c in df.columns
            ]
            fh.write(df[cols].to_markdown(index=False))
            fh.write("\n\n")
    console.print(f"[green]✓ {full_csv}[/green]")
    console.print(f"[green]✓ {md_path}[/green]")
    if not splits["conserved"].empty:
        top = ", ".join(splits["conserved"]["tf"].head(10).tolist())
        console.print(f"[bold green]Conserved top-{cap} (head): {top} …[/bold green]")


@cli.command("score-baseline")
@click.option(
    "--resting-celltypes",
    default="stromal_fibroblast",
    show_default=True,
    help="Comma-separated cell_type labels pooled as the resting baseline.",
)
@click.option(
    "--decidualized-celltypes",
    default="pre_decidual_stromal,decidual_stromal,senescent_decidual",
    show_default=True,
    help="Comma-separated cell_type labels pooled as the decidualized end-state.",
)
@click.option("--min-cells", default=20, show_default=True, type=int)
def score_baseline_cmd(
    resting_celltypes: str, decidualized_celltypes: str, min_cells: int
) -> None:
    """Baseline-priming test (Q4.1): per-species priming distance from
    resting → decidualized stroma, plus between-species comparison at
    the resting baseline.
    """
    from pathlib import Path

    import anndata as ad

    from src.scoring.baseline_priming import (
        BaselinePrimingConfig,
        baseline_priming,
        between_species_resting,
    )

    project_root = Path(__file__).resolve().parent.parent
    scored_path = project_root / "results" / "scored" / "stromal_scored.h5ad"
    if not scored_path.exists():
        console.print(f"[red]Missing {scored_path}. Run `wombat score-decidua` first.[/red]")
        raise SystemExit(1)

    console.print(f"[blue]Loading {scored_path}...[/blue]")
    adata = ad.read_h5ad(scored_path)
    score_cols = [c for c in adata.obs.columns if c.endswith("_score")]
    if not score_cols:
        console.print("[red]No *_score columns found in adata.obs.[/red]")
        raise SystemExit(1)

    cfg = BaselinePrimingConfig(
        resting_celltypes=tuple(s.strip() for s in resting_celltypes.split(",") if s.strip()),
        decidualized_celltypes=tuple(
            s.strip() for s in decidualized_celltypes.split(",") if s.strip()
        ),
        min_cells_per_group=min_cells,
    )
    console.print(
        f"[blue]Scoring {len(score_cols)} modules: "
        f"resting={cfg.resting_celltypes}, decidualized={cfg.decidualized_celltypes}...[/blue]"
    )
    priming = baseline_priming(adata, score_cols=score_cols, config=cfg)
    between = between_species_resting(adata, score_cols=score_cols, config=cfg)

    out_dir = project_root / "results" / "reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    priming_csv = out_dir / "baseline_priming.csv"
    between_csv = out_dir / "baseline_priming_between_species.csv"
    md_path = out_dir / "baseline_priming.md"
    priming.to_csv(priming_csv, index=False)
    between.to_csv(between_csv, index=False)

    # Hypothesis-1 verdict for the decidual_score row(s).
    verdict_lines: list[str] = []
    if "decidual_score" in priming["score"].values:
        ds = priming[priming["score"] == "decidual_score"].set_index("species")
        if {"human", "mouse"}.issubset(ds.index):
            h_d = ds.loc["human", "priming_distance"]
            m_d = ds.loc["mouse", "priming_distance"]
            gap = m_d - h_d
            verdict_lines.append(
                f"- Human priming distance: **{h_d:.3f}** (n_resting="
                f"{int(ds.loc['human', 'n_resting'])}, "
                f"n_decidualized={int(ds.loc['human', 'n_decidualized'])})"
            )
            verdict_lines.append(
                f"- Mouse priming distance: **{m_d:.3f}** (n_resting="
                f"{int(ds.loc['mouse', 'n_resting'])}, "
                f"n_decidualized={int(ds.loc['mouse', 'n_decidualized'])})"
            )
            verdict_lines.append(f"- Gap (mouse − human): **{gap:+.3f}**")
            if gap > 0.5:
                verdict_lines.append(
                    "- **Supports hypothesis 1** (lowered activation threshold): "
                    "human resting stroma is closer to the decidualized end-state."
                )
            elif gap < -0.5:
                verdict_lines.append(
                    "- **Refutes hypothesis 1** in the expected direction: mouse "
                    "resting stroma is closer to the end-state than human."
                )
            else:
                verdict_lines.append(
                    "- **Inconclusive / refutes hypothesis 1**: the two species "
                    "have comparable priming distances (|gap| ≤ 0.5 sd). "
                    "Focus Q4.2 / Q4.3 on hypotheses 2 and 4."
                )

    with open(md_path, "w") as fh:
        fh.write("# Baseline-priming test (Q4.1)\n\n")
        fh.write(
            "Hypothesis 1 of the Q4 convergent-evolution question: spontaneous-\n"
            "deciduator stroma sits at a **lowered activation threshold**. If true,\n"
            "human (spontaneous) resting stromal cells should already score higher\n"
            f"on `decidual_score` than mouse (induced) resting cells, and the\n"
            "within-species priming distance (Cohen's d, resting → decidualized)\n"
            "should be *smaller* in the spontaneous species.\n\n"
            f"- Resting celltype(s): `{', '.join(cfg.resting_celltypes)}`\n"
            f"- Decidualized celltype(s): `{', '.join(cfg.decidualized_celltypes)}`\n"
            f"- Min cells per group: {min_cells}\n\n"
        )
        if verdict_lines:
            fh.write("## Hypothesis 1 verdict (decidual_score)\n\n")
            fh.write("\n".join(verdict_lines))
            fh.write("\n\n")
        fh.write("## Per-species priming distance (all modules)\n\n")
        fh.write(priming.to_markdown(index=False, floatfmt=".3f"))
        fh.write("\n\n## Between-species resting comparison\n\n")
        fh.write(
            "`cohens_d_b_minus_a` > 0 means species_b has higher resting score "
            "than species_a. For `decidual_score` with species_a=human, "
            "species_b=mouse, a **negative** value supports hypothesis 1.\n\n"
        )
        fh.write(between.to_markdown(index=False, floatfmt=".3f"))
        fh.write("\n")

    console.print(f"[green]✓ {priming_csv}[/green]")
    console.print(f"[green]✓ {between_csv}[/green]")
    console.print(f"[green]✓ {md_path}[/green]")
    for line in verdict_lines:
        console.print(line)


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
