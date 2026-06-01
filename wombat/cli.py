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
@click.option(
    "--target",
    default="mouse",
    show_default=True,
    help="Target species name (must exist in configs/species.yaml). "
    "Ignored when --all-tier-b is set.",
)
@click.option(
    "--source",
    default="human",
    show_default=True,
    help="Source species name (must exist in configs/species.yaml).",
)
@click.option(
    "--all-tier-b",
    is_flag=True,
    help="Build a backbone for every Tier B species (continues past per-target failures).",
)
def orthologs_build(no_gprofiler: bool, target: str, source: str, all_tier_b: bool) -> None:
    """Build human→target ortholog backbone(s).

    Default: human → mouse, written to
    ``results/orthologs/backbone.parquet`` (preserves the historical
    output path used by downstream consumers).

    With ``--target <name>`` the backbone is written to
    ``results/orthologs/backbone__<source>_<target>.parquet``.

    With ``--all-tier-b`` every species in ``configs/species.yaml`` with
    ``tier: B`` becomes a target in turn; per-target failures are logged
    and the loop continues so a single missing Ensembl dataset (e.g.
    tenrec) does not block the rest.
    """
    from pathlib import Path

    from src.orthologs.backbone import build_backbone
    from wombat.config import load_config

    project_root = Path(__file__).resolve().parent.parent
    cache_dir = project_root / "results" / "orthologs" / "cache"
    orth_dir = project_root / "results" / "orthologs"

    if all_tier_b:
        targets = [s["name"] for s in load_config("species") if s.get("tier") == "B"]
        if not targets:
            console.print("[red]No Tier B species found in configs/species.yaml[/red]")
            raise SystemExit(1)
    else:
        targets = [target]

    successes: list[tuple[str, int]] = []
    failures: list[tuple[str, str]] = []

    for tgt in targets:
        if not all_tier_b and tgt == "mouse" and source == "human":
            output_path = orth_dir / "backbone.parquet"  # historical default
        else:
            output_path = orth_dir / f"backbone__{source}_{tgt}.parquet"

        console.print(f"[blue]Building ortholog backbone ({source} → {tgt})...[/blue]")
        try:
            backbone = build_backbone(
                source=source,
                target=tgt,
                cache_dir=cache_dir,
                output_path=output_path,
                use_gprofiler=not no_gprofiler,
            )
        except Exception as exc:  # noqa: BLE001 — per-target continue
            console.print(f"[red]✗ {source} → {tgt} failed: {exc}[/red]")
            failures.append((tgt, str(exc)))
            continue

        successes.append((tgt, len(backbone)))
        console.print(f"[green]✓ {source} → {tgt}: {len(backbone)} rows → {output_path}[/green]")

    if all_tier_b:
        console.print(f"[bold]Summary: {len(successes)} succeeded, {len(failures)} failed[/bold]")
        for tgt, n in successes:
            console.print(f"  ✓ {tgt}: {n} rows")
        for tgt, err in failures:
            console.print(f"  ✗ {tgt}: {err}")


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


@cli.command("trait-contrast")
@click.option("--min-module-genes", default=3, show_default=True, type=int)
@click.option("--min-samples-per-arm", default=2, show_default=True, type=int)
def trait_contrast_cmd(min_module_genes: int, min_samples_per_arm: int) -> None:
    """Trait contrast (Q4.2): pseudobulk module amplitude in spontaneous
    vs induced deciduators, across the GSE155170 gene-level bulk deposits.

    CAVEAT: in this four-species subset the trait is perfectly confounded
    with clade (catarrhine/bat vs rodent), so the contrast measures a
    trait-or-clade difference, not a phylogeny-controlled trait effect.
    """
    from pathlib import Path

    import anndata as ad
    import pyarrow.parquet as pq

    from src.scoring.trait_contrast import score_species_pseudobulk, trait_contrast
    from wombat.config import load_config

    project_root = Path(__file__).resolve().parent.parent
    processed_dir = project_root / "results" / "processed"
    orth_dir = project_root / "results" / "orthologs"

    gene_h5ads = sorted(processed_dir.glob("*__gene.h5ad"))
    if not gene_h5ads:
        console.print(f"[red]No *__gene.h5ad files in {processed_dir}.[/red]")
        raise SystemExit(1)

    markers_cfg = load_config("markers")
    gene_sets = markers_cfg.get("score_gene_sets") if isinstance(markers_cfg, dict) else None
    if not gene_sets:
        console.print("[red]No score_gene_sets in configs/markers.yaml.[/red]")
        raise SystemExit(1)

    species_cfg = {s["name"]: s for s in load_config("species")}

    per_sample_frames: list = []
    coverage_rows: list[dict[str, object]] = []
    trait_by_species: dict[str, bool] = {}
    clade_by_species: dict[str, str] = {}

    for h5ad_path in gene_h5ads:
        # Filename convention: <DATASET>__<species>__gene.h5ad
        parts = h5ad_path.name.split("__")
        species = parts[1] if len(parts) >= 3 else h5ad_path.stem
        meta = species_cfg.get(species)
        if meta is None or "spontaneous_decidualization" not in meta:
            console.print(f"[yellow]  Skipping {species}: no trait label in species.yaml[/yellow]")
            continue
        trait_by_species[species] = bool(meta["spontaneous_decidualization"])
        clade_by_species[species] = str(meta.get("clade", "?"))

        per_species = orth_dir / f"backbone__human_{species}.parquet"
        if per_species.exists():
            backbone_df = pq.read_table(per_species).to_pandas()
        else:
            # Tier-A species (e.g. mouse) live in the main backbone.
            full = pq.read_table(orth_dir / "backbone.parquet").to_pandas()
            backbone_df = full[full["target_species"] == species].reset_index(drop=True)
        if not len(backbone_df):
            console.print(
                f"[yellow]  No backbone rows for {species}; symbol-only mapping.[/yellow]"
            )
            backbone_df = None

        console.print(f"[blue]Scoring {species} ({h5ad_path.name})...[/blue]")
        adata = ad.read_h5ad(h5ad_path)
        scored = score_species_pseudobulk(
            adata,
            gene_sets=gene_sets,
            species=species,
            backbone_df=backbone_df,
            min_module_genes=min_module_genes,
        )
        per_sample_frames.append(scored)

        cov = scored.drop_duplicates("score")[["score", "n_mapped", "n_by_id", "n_by_symbol"]]
        for _, r in cov.iterrows():
            coverage_rows.append(
                {
                    "species": species,
                    "clade": clade_by_species[species],
                    "trait_positive": trait_by_species[species],
                    "score": r["score"],
                    "n_mapped": int(r["n_mapped"]),
                    "n_by_id": int(r["n_by_id"]),
                    "n_by_symbol": int(r["n_by_symbol"]),
                }
            )

    if not per_sample_frames:
        console.print("[red]No species scored.[/red]")
        raise SystemExit(1)

    import pandas as pd

    scores = pd.concat(per_sample_frames, ignore_index=True)
    coverage = pd.DataFrame(coverage_rows)
    contrast = trait_contrast(
        scores, trait_by_species=trait_by_species, min_samples_per_arm=min_samples_per_arm
    )

    out_dir = project_root / "results" / "reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    scores_csv = out_dir / "trait_contrast_scores.csv"
    coverage_csv = out_dir / "trait_contrast_coverage.csv"
    contrast_csv = out_dir / "trait_contrast.csv"
    md_path = out_dir / "trait_contrast.md"
    scores.to_csv(scores_csv, index=False)
    coverage.to_csv(coverage_csv, index=False)
    contrast.to_csv(contrast_csv, index=False)

    pos = sorted(s for s, t in trait_by_species.items() if t)
    neg = sorted(s for s, t in trait_by_species.items() if not t)
    pos_clades = sorted({clade_by_species[s] for s in pos})
    neg_clades = sorted({clade_by_species[s] for s in neg})
    confounded = not (set(pos_clades) & set(neg_clades))

    with open(md_path, "w") as fh:
        fh.write("# Trait contrast (Q4.2)\n\n")
        fh.write(
            "Does the conserved decidual program reach a **higher pseudobulk\n"
            "amplitude** in species that decidualize spontaneously than in those\n"
            "that require embryonic induction? Each module score is the mean\n"
            "within-sample z-score (CPM-log1p) of its mapped genes; modules are\n"
            "mapped from human symbols through the ortholog backbone with a\n"
            "gene-symbol fallback (recovers markers whose gene IDs drifted between\n"
            "annotation releases).\n\n"
        )
        fh.write(
            f"- Trait-positive (spontaneous): "
            f"{', '.join(f'{s} [{clade_by_species[s]}]' for s in pos)}\n"
        )
        fh.write(
            f"- Trait-negative (induced): "
            f"{', '.join(f'{s} [{clade_by_species[s]}]' for s in neg)}\n\n"
        )
        if confounded:
            fh.write(
                "> **Confound caveat.** In this subset the trait is perfectly\n"
                f"> confounded with clade (positive={pos_clades}, "
                f"negative={neg_clades}); the contrast cannot separate a trait\n"
                "> effect from a clade effect. A phylogeny-controlled test needs a\n"
                "> within-clade trait contrast (e.g. a trait-negative catarrhine or\n"
                "> a trait-positive rodent) that the current data do not provide.\n\n"
            )
        fh.write("## Trait contrast (spontaneous − induced)\n\n")
        fh.write(
            "`delta_pos_minus_neg` > 0 and `fdr` < 0.05 ⇒ the module is more\n"
            "strongly expressed in spontaneous deciduators.\n\n"
        )
        fh.write(contrast.to_markdown(index=False, floatfmt=".3f"))
        fh.write("\n\n## Per-species module coverage\n\n")
        fh.write(
            "`n_by_symbol` > 0 marks genes recovered through the symbol fallback "
            "after the gene-ID join missed them.\n\n"
        )
        cov_wide = coverage.pivot_table(
            index="score", columns="species", values="n_mapped", aggfunc="first"
        )
        fh.write(cov_wide.to_markdown(floatfmt=".0f"))
        fh.write("\n")

    console.print(f"[green]✓ {scores_csv}[/green]")
    console.print(f"[green]✓ {coverage_csv}[/green]")
    console.print(f"[green]✓ {contrast_csv}[/green]")
    console.print(f"[green]✓ {md_path}[/green]")
    sig = contrast[(contrast["fdr"] < 0.05) & (contrast["delta_pos_minus_neg"] > 0)]
    if len(sig):
        console.print(
            f"[green]{len(sig)} module(s) higher in spontaneous deciduators (FDR<0.05):[/green]"
        )
        for _, r in sig.iterrows():
            console.print(
                f"  - {r['score']}: Δ={r['delta_pos_minus_neg']:+.3f}, d={r['cohens_d']:+.2f}, FDR={r['fdr']:.3f}"
            )
    else:
        console.print("[yellow]No module passed FDR<0.05 in the spontaneous direction.[/yellow]")


@cli.command("cis-regulatory")
@click.option(
    "--window",
    default=50_000,
    show_default=True,
    type=int,
    help="±bp around each decidual-gene TSS for the proximity test.",
)
def cis_regulatory_cmd(window: int) -> None:
    """Cis-regulatory TE overlap (Q4.3): test the Lynch/Wagner hypothesis
    that the human decidual enhancer landscape is transposable-element
    derived, using processed GSE61793 ChIP/DNase peaks (hg19) and UCSC
    RepeatMasker.

    Two readouts: (1) a genome-wide TE census of all peaks per assay, and
    (2) a decidual-gene proximity test asking whether the flagged ancient
    families (MER20/MER41) are enriched in enhancers near the decidual
    module genes. CAVEAT: GSE61793 is human-only ChIP, so this is a
    descriptive landscape test, not a cross-species trait contrast.
    """
    from pathlib import Path

    import pandas as pd

    from src.cis_regulatory.genes import gene_windows, load_tss, matched_symbols
    from src.cis_regulatory.peaks import load_all_peaks, peak_qc
    from src.cis_regulatory.te_overlap import (
        load_rmsk,
        near_gene_te_enrichment,
        te_enrichment_all,
    )
    from wombat.config import load_config

    project_root = Path(__file__).resolve().parent.parent
    peaks_dir = project_root / "results" / "raw" / "GSE61793"
    ref_dir = project_root / "results" / "raw" / "reference"
    rmsk_path = ref_dir / "rmsk_hg19.txt.gz"
    refgene_path = ref_dir / "refGene_hg19.txt.gz"

    peaks = load_all_peaks(peaks_dir)
    if not peaks:
        console.print(f"[red]No GSE61793 peak BEDs in {peaks_dir}.[/red]")
        raise SystemExit(1)
    for path in (rmsk_path, refgene_path):
        if not path.exists():
            console.print(f"[red]Missing reference: {path}[/red]")
            raise SystemExit(1)

    qc = pd.DataFrame(peak_qc(df, assay) for assay, df in peaks.items())

    console.print("[blue]Loading RepeatMasker (hg19)...[/blue]")
    rmsk = load_rmsk(rmsk_path)
    console.print("[blue]Computing genome-wide TE census...[/blue]")
    census = te_enrichment_all(peaks, rmsk)

    markers_cfg = load_config("markers")
    gene_sets = markers_cfg.get("score_gene_sets") if isinstance(markers_cfg, dict) else None
    if not gene_sets:
        console.print("[red]No score_gene_sets in configs/markers.yaml.[/red]")
        raise SystemExit(1)
    genes = sorted({g for v in gene_sets.values() for g in v})

    tss = load_tss(refgene_path)
    matched = matched_symbols(tss, genes)
    windows = gene_windows(tss, genes, window=window)
    console.print(
        f"[blue]Decidual-gene proximity: {len(matched)}/{len(genes)} symbols, "
        f"{len(windows)} merged windows (±{window:,} bp).[/blue]"
    )
    proximity = pd.concat(
        [near_gene_te_enrichment(df, rmsk, windows, assay) for assay, df in peaks.items()],
        ignore_index=True,
    )

    out_dir = project_root / "results" / "reports"
    out_dir.mkdir(parents=True, exist_ok=True)
    qc_csv = out_dir / "cis_regulatory_peak_qc.csv"
    census_csv = out_dir / "cis_regulatory_te_census.csv"
    proximity_csv = out_dir / "cis_regulatory_proximity.csv"
    md_path = out_dir / "cis_regulatory.md"
    qc.to_csv(qc_csv, index=False)
    census.to_csv(census_csv, index=False)
    proximity.to_csv(proximity_csv, index=False)

    sig = proximity[proximity["fisher_p_greater"] < 0.05]
    with open(md_path, "w") as fh:
        fh.write("# Cis-regulatory TE overlap (Q4.3)\n\n")
        fh.write(
            "Is the human decidual *cis*-regulatory landscape disproportionately\n"
            "**transposable-element derived**, as the Lynch/Wagner model of\n"
            "co-opted endometrial regulation predicts? We overlap processed\n"
            "GSE61793 peaks (hg19; H3K27ac active enhancers, H3K4me3 promoters,\n"
            "DNaseI open chromatin) with UCSC RepeatMasker and flag the ancient\n"
            "families MER20 and MER41 implicated in that model.\n\n"
        )
        fh.write(
            "> **Scope caveat.** GSE61793 is human-only ChIP, so this is a\n"
            "> descriptive landscape test, **not** a cross-species trait contrast.\n"
            "> The spontaneous-vs-induced comparison the hypothesis ultimately\n"
            "> needs is not answerable from a single species' regulatory maps.\n\n"
        )
        fh.write("## Peak sets\n\n")
        fh.write(qc.to_markdown(index=False, floatfmt=".0f"))
        fh.write("\n\n## Genome-wide TE census\n\n")
        fh.write(
            "Fraction of peaks (per assay) overlapping each RepeatMasker\n"
            "category. SVA retroposons are filed under class `Other` in hg19.\n\n"
        )
        fh.write(census.to_markdown(index=False, floatfmt=".4f"))
        fh.write("\n\n## Decidual-gene proximity test (the Lynch prediction)\n\n")
        fh.write(
            f"Among peaks within ±{window:,} bp of a decidual module gene TSS\n"
            f"({len(matched)}/{len(genes)} symbols, {len(windows)} merged windows),\n"
            "is a flagged family enriched vs peaks elsewhere? One-sided Fisher's\n"
            'exact test (`alternative="greater"`).\n\n'
        )
        fh.write(proximity.to_markdown(index=False, floatfmt=".4f"))
        fh.write("\n\n## Verdict\n\n")
        if len(sig):
            fh.write(
                "Flagged TE families **enriched** in decidual-gene peaks (Fisher p<0.05):\n\n"
            )
            for _, r in sig.iterrows():
                fh.write(
                    f"- {r['assay']} / {r['family']}: "
                    f"{r['near_hit']}/{r['n_near']} near vs "
                    f"{r['far_hit']}/{r['n_far']} far "
                    f"(OR={r['odds_ratio']:.2f}, p={r['fisher_p_greater']:.3g})\n"
                )
        else:
            fh.write(
                "No flagged family (MER20/MER41) is significantly enriched in\n"
                "decidual-gene peaks at p<0.05. With the project's compact\n"
                "decidual panel the proximity test is **underpowered** (few\n"
                "near-gene peaks vs a <1% family base rate); the genome-wide\n"
                "census still shows the bulk of enhancers are TE-derived, but\n"
                "the gene-specific Lynch signal is not recovered from this\n"
                "human-only ChIP.\n"
            )
        fh.write("\n")

    console.print(f"[green]✓ {qc_csv}[/green]")
    console.print(f"[green]✓ {census_csv}[/green]")
    console.print(f"[green]✓ {proximity_csv}[/green]")
    console.print(f"[green]✓ {md_path}[/green]")
    if len(sig):
        console.print(
            f"[green]{len(sig)} flagged family enrichment(s) in decidual-gene peaks "
            "(Fisher p<0.05).[/green]"
        )
    else:
        console.print(
            "[yellow]No MER20/MER41 enrichment in decidual-gene peaks "
            "(underpowered panel).[/yellow]"
        )


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
