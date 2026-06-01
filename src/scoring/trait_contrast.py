"""Multi-species trait contrast (Q4.2).

Q4.2 asks whether the conserved decidual program is *more strongly
expressed* in species that decidualize **spontaneously** (menstruating
catarrhines, phyllostomid bats) than in species that only decidualize on
embryonic induction (most rodents). The Q4.1 within-atlas
resting→decidualized priming distance does not transfer to these bulk
endometrial deposits (different tissue, no single-cell resting/decidual
split), so the Q4.2 readout is a **pseudobulk module-amplitude contrast**:

1. Each species ships a gene-level bulk AnnData (one row per sample) with
   ``ENS<spp>G`` gene IDs in ``var_names`` and a ``gene_symbol`` column
   (produced by :mod:`src.ingest.tx2gene`).
2. The 8 conserved scoring modules (``configs/markers.yaml``) are defined
   on human symbols. For each species they are mapped to the species'
   own gene space through the ortholog backbone, with a **gene-symbol
   fallback** — essential because some deposits (notably baboon
   GSE155170) were annotated against an older Ensembl release whose gene
   IDs have since drifted, so an ID-only join silently drops canonical
   markers (IGFBP1, CCND1, ...) that are still recoverable by symbol.
3. Each sample gets a within-sample, z-scored module score (mean z of the
   mapped module genes), making scores comparable across samples,
   libraries, and species.
4. Module scores are contrasted between the trait-positive and
   trait-negative species pools (Welch's t, Cohen's d, Benjamini-Hochberg
   FDR across modules).

The pure functions here take already-loaded AnnData / backbone DataFrames
so they are unit-testable without IO; the ``wombat trait-contrast`` CLI
command does the file loading and report writing.

**Caveat carried from the scoping work:** in the GSE155170 four-species
subset the trait (spontaneous vs induced) is perfectly confounded with
clade (catarrhine/bat vs rodent), so the contrast measures a
trait-or-clade difference, not a phylogeny-controlled trait effect. A
true phylogenetic regression needs a within-clade trait contrast that the
current data do not provide; this is documented, not corrected here.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from anndata import AnnData
from scipy import stats

from src.scoring.baseline_priming import _cohens_d
from src.scoring.null import _bh


@dataclass(frozen=True)
class GeneMapping:
    """Result of mapping one human module to a species' gene space.

    ``by_gene`` keeps the per-human-gene grouping (one human symbol may
    map to several paralogous target genes under one-to-many orthology)
    so the scorer can collapse paralogs to a single representative and
    avoid letting paralog count bias the module amplitude.
    """

    by_gene: dict[str, list[str]]
    n_by_id: int
    n_by_symbol: int

    @property
    def matched_var_names(self) -> list[str]:
        """Flat, de-duplicated list of all matched target var_names."""
        out: list[str] = []
        seen: set[str] = set()
        for vns in self.by_gene.values():
            for vn in vns:
                if vn not in seen:
                    seen.add(vn)
                    out.append(vn)
        return out


def map_module_genes(
    human_symbols: list[str],
    var_names,
    var_symbols: dict[str, str] | None = None,
    backbone_df: pd.DataFrame | None = None,
) -> GeneMapping:
    """Map human gene symbols onto a species' ``var_names``.

    Two stages, applied per human symbol:

    1. **Gene-ID join** — via the ortholog backbone's
       ``source_symbol → target_gene_id`` (Tier-1 rows preferred), kept
       only when the target gene ID is present in ``var_names``.
    2. **Symbol fallback** — when the ID join finds nothing, match the
       human symbol (and the backbone's ``target_symbol`` for it) against
       the species' ``gene_symbol`` column (case-insensitive). This
       recovers genes whose Ensembl *gene ID* drifted between the
       deposit's annotation release and the backbone's release but whose
       symbol is stable.

    Parameters
    ----------
    human_symbols : list[str]
        Human gene symbols defining the module.
    var_names : iterable
        The species AnnData's ``var_names`` (gene IDs, or human symbols
        when ``backbone_df`` is None).
    var_symbols : dict[str, str], optional
        Mapping ``var_name → gene_symbol`` for the symbol fallback.
    backbone_df : pd.DataFrame, optional
        Ortholog backbone filtered to the target species (columns
        ``source_symbol``, ``target_gene_id``, ``target_symbol``,
        ``tier``). ``None`` means the AnnData is already in human-symbol
        space (e.g. human itself) and symbols are matched directly.

    Returns
    -------
    GeneMapping
        Per-human-gene matched ``var_names`` plus how many human genes
        were resolved by ID vs by symbol (for coverage reporting).
    """
    var_list = [str(v) for v in var_names]
    var_set = set(var_list)
    sym_to_var: dict[str, list[str]] = {}
    if var_symbols:
        for vn, sym in var_symbols.items():
            if sym:
                sym_to_var.setdefault(str(sym).upper(), []).append(str(vn))

    # Human symbol → backbone target gene IDs / target symbols.
    id_index: dict[str, list[str]] = {}
    tgtsym_index: dict[str, list[str]] = {}
    if backbone_df is not None and len(backbone_df):
        df = backbone_df
        if "tier" in df.columns:
            df = df.sort_values("tier")  # Tier-1 first so it wins ties
        for src, tid, tsym in zip(
            df["source_symbol"],
            df["target_gene_id"],
            df.get("target_symbol", df["target_gene_id"]),
            strict=False,
        ):
            if not src:
                continue
            su = str(src).upper()
            if tid:
                id_index.setdefault(su, []).append(str(tid))
            if isinstance(tsym, str) and tsym:
                tgtsym_index.setdefault(su, []).append(tsym.upper())

    by_gene: dict[str, list[str]] = {}
    n_by_id = 0
    n_by_symbol = 0
    for g in human_symbols:
        gu = g.upper()
        # Stage 1: gene-ID join.
        id_hits = list(dict.fromkeys(tid for tid in id_index.get(gu, []) if tid in var_set))
        if id_hits:
            by_gene[g] = id_hits
            n_by_id += 1
            continue
        # Stage 2: symbol fallback (human symbol + backbone target symbols).
        candidate_syms = {gu, *tgtsym_index.get(gu, [])}
        sym_hits: list[str] = []
        if sym_to_var:
            for cs in candidate_syms:
                sym_hits.extend(sym_to_var.get(cs, []))
        elif backbone_df is None:
            # Human-symbol-indexed AnnData: var_names themselves are symbols.
            sym_hits = [v for v in var_list if v.upper() in candidate_syms]
        sym_hits = list(dict.fromkeys(sym_hits))
        if sym_hits:
            by_gene[g] = sym_hits
            n_by_symbol += 1

    return GeneMapping(by_gene=by_gene, n_by_id=n_by_id, n_by_symbol=n_by_symbol)


def _within_sample_z(adata: AnnData) -> np.ndarray:
    """CPM-log1p normalise then z-score each sample across genes.

    Returns a dense ``(n_samples, n_genes)`` array whose rows are each
    sample's per-gene z-scores, so a module score (mean z over the
    module's genes) is comparable across samples and species regardless
    of library size or absolute expression scale.
    """
    x = adata.X
    x = x.toarray() if hasattr(x, "toarray") else np.asarray(x)
    x = np.asarray(x, dtype=float)
    libsize = x.sum(axis=1, keepdims=True)
    libsize[libsize == 0] = 1.0
    cpm = x / libsize * 1e6
    logx = np.log1p(cpm)
    mean = logx.mean(axis=1, keepdims=True)
    std = logx.std(axis=1, ddof=0, keepdims=True)
    std[std == 0] = 1.0
    return (logx - mean) / std


def score_species_pseudobulk(
    adata: AnnData,
    gene_sets: dict[str, list[str]],
    species: str,
    backbone_df: pd.DataFrame | None = None,
    min_module_genes: int = 3,
) -> pd.DataFrame:
    """Per-sample module scores for one species' bulk AnnData.

    Each module score is the mean within-sample z-score of the module's
    mapped genes, where one-to-many orthologs are collapsed to a single
    representative per human gene (the most-expressed paralog in that
    sample) so paralog count does not bias the amplitude. ``n_mapped`` is
    therefore the number of **human genes** resolved. Modules mapping
    fewer than ``min_module_genes`` human genes are emitted with NaN
    scores (and the true ``n_mapped``) so the caller can see — and
    exclude — under-covered modules rather than trusting a one- or
    two-gene proxy.

    Returns one row per ``(sample, score)`` with columns ``species``,
    ``sample``, ``score``, ``value``, ``n_mapped``, ``n_by_id``,
    ``n_by_symbol``.
    """
    var_symbols = None
    if "gene_symbol" in adata.var.columns:
        var_symbols = dict(
            zip(
                adata.var_names.astype(str),
                adata.var["gene_symbol"].astype(str),
                strict=False,
            )
        )

    z = _within_sample_z(adata)
    var_pos = {str(v): i for i, v in enumerate(adata.var_names)}
    sample_ids = [str(s) for s in adata.obs_names]

    rows: list[dict[str, object]] = []
    for score_name, genes in gene_sets.items():
        mapping = map_module_genes(genes, adata.var_names, var_symbols, backbone_df)
        # Collapse paralogs: one representative z per human gene = the
        # most-expressed paralog (max z) in each sample.
        per_gene_z: list[np.ndarray] = []
        for var_list in mapping.by_gene.values():
            idx = [var_pos[v] for v in var_list if v in var_pos]
            if idx:
                per_gene_z.append(z[:, idx].max(axis=1))
        n_mapped = len(per_gene_z)
        enough = n_mapped >= min_module_genes
        if enough and per_gene_z:
            module_z = np.vstack(per_gene_z).mean(axis=0)
        else:
            module_z = np.full(z.shape[0], np.nan)
        for sample, val in zip(sample_ids, module_z, strict=False):
            rows.append(
                {
                    "species": species,
                    "sample": sample,
                    "score": score_name,
                    "value": float(val),
                    "n_mapped": n_mapped,
                    "n_by_id": mapping.n_by_id,
                    "n_by_symbol": mapping.n_by_symbol,
                }
            )
    return pd.DataFrame(rows)


CONTRAST_COLS = (
    "score",
    "n_species_pos",
    "n_species_neg",
    "n_samples_pos",
    "n_samples_neg",
    "mean_pos",
    "mean_neg",
    "delta_pos_minus_neg",
    "cohens_d",
    "welch_t",
    "welch_p",
    "fdr",
)


def trait_contrast(
    scores: pd.DataFrame,
    trait_by_species: dict[str, bool],
    min_samples_per_arm: int = 2,
) -> pd.DataFrame:
    """Contrast per-module scores between trait-positive and -negative pools.

    Parameters
    ----------
    scores : pd.DataFrame
        Long per-sample scores (output of
        :func:`score_species_pseudobulk`, concatenated across species).
    trait_by_species : dict[str, bool]
        ``species → spontaneous_decidualization`` flag.
    min_samples_per_arm : int
        Minimum samples in each arm for a module to be tested; otherwise
        the statistics are NaN.

    Returns
    -------
    pd.DataFrame
        One row per module: arm sizes, arm means, ``delta`` (pos − neg),
        Cohen's d (pos vs neg), Welch's t / p, and Benjamini-Hochberg
        ``fdr`` across the tested modules. ``delta > 0`` means the
        spontaneous-deciduator pool expresses the module more strongly.
    """
    required = {"species", "sample", "score", "value"}
    missing = required - set(scores.columns)
    if missing:
        msg = f"scores missing required columns: {sorted(missing)}"
        raise ValueError(msg)

    work = scores.copy()
    work["trait_pos"] = work["species"].map(trait_by_species)
    if work["trait_pos"].isna().any():
        unknown = sorted(work.loc[work["trait_pos"].isna(), "species"].unique())
        msg = f"trait_by_species has no entry for species: {unknown}"
        raise ValueError(msg)

    rows: list[dict[str, object]] = []
    for score_name in work["score"].drop_duplicates():
        sub = work[work["score"] == score_name]
        pos = sub[sub["trait_pos"]]
        neg = sub[~sub["trait_pos"]]
        pos_vals = pos["value"].to_numpy(dtype=float)
        neg_vals = neg["value"].to_numpy(dtype=float)
        pos_vals = pos_vals[~np.isnan(pos_vals)]
        neg_vals = neg_vals[~np.isnan(neg_vals)]

        if len(pos_vals) < min_samples_per_arm or len(neg_vals) < min_samples_per_arm:
            t_stat, p_val, d = float("nan"), float("nan"), float("nan")
        else:
            t_stat, p_val = stats.ttest_ind(pos_vals, neg_vals, equal_var=False)
            d = _cohens_d(neg_vals, pos_vals)  # d > 0 ⇒ pos mean > neg mean

        mean_pos = float(np.mean(pos_vals)) if len(pos_vals) else float("nan")
        mean_neg = float(np.mean(neg_vals)) if len(neg_vals) else float("nan")
        rows.append(
            {
                "score": score_name,
                "n_species_pos": int(pos["species"].nunique()),
                "n_species_neg": int(neg["species"].nunique()),
                "n_samples_pos": int(len(pos_vals)),
                "n_samples_neg": int(len(neg_vals)),
                "mean_pos": mean_pos,
                "mean_neg": mean_neg,
                "delta_pos_minus_neg": mean_pos - mean_neg,
                "cohens_d": d,
                "welch_t": float(t_stat),
                "welch_p": float(p_val),
            }
        )

    out = pd.DataFrame(rows, columns=[c for c in CONTRAST_COLS if c != "fdr"])
    pvals = out["welch_p"].to_numpy(dtype=float)
    fdr = np.full(len(pvals), np.nan)
    tested = ~np.isnan(pvals)
    if tested.any():
        fdr[tested] = _bh(pvals[tested])
    out["fdr"] = fdr
    return out
