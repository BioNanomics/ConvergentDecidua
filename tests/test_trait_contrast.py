"""Unit tests for src/scoring/trait_contrast.py (Q4.2)."""

from __future__ import annotations

import numpy as np
import pandas as pd
from anndata import AnnData

from src.scoring.trait_contrast import (
    map_module_genes,
    score_species_pseudobulk,
    trait_contrast,
)


def _backbone() -> pd.DataFrame:
    """Tiny ortholog backbone (human → some species)."""
    return pd.DataFrame(
        {
            "source_symbol": ["GENEA", "GENEB", "GENEC", "GENED"],
            # GENEA's ID is present in the matrix (ID join).
            # GENEB's backbone ID drifted away (recover by human symbol).
            # GENEC's backbone ID drifted; recover via backbone target_symbol.
            # GENED maps to nothing.
            "target_gene_id": ["ENSX1", "ENSX_OLD_B", "ENSX_OLD_C", "ENSX_GONE"],
            "target_symbol": ["Asym", "Bsym", "Csym", "Dsym"],
            "tier": [1, 1, 1, 1],
        }
    )


def test_map_module_genes_id_join():
    bb = _backbone()
    var_names = ["ENSX1", "ENSX2", "ENSX3"]
    mapping = map_module_genes(["GENEA"], var_names, var_symbols=None, backbone_df=bb)
    assert mapping.matched_var_names == ["ENSX1"]
    assert mapping.n_by_id == 1
    assert mapping.n_by_symbol == 0


def test_map_module_genes_symbol_fallback_human_symbol():
    bb = _backbone()
    var_names = ["ENSX1", "ENSX2", "ENSX3"]
    # ENSX2 carries the human symbol GENEB even though the backbone ID drifted.
    var_symbols = {"ENSX1": "Asym", "ENSX2": "GENEB", "ENSX3": "other"}
    mapping = map_module_genes(["GENEB"], var_names, var_symbols=var_symbols, backbone_df=bb)
    assert mapping.matched_var_names == ["ENSX2"]
    assert mapping.n_by_id == 0
    assert mapping.n_by_symbol == 1


def test_map_module_genes_symbol_fallback_backbone_target_symbol():
    bb = _backbone()
    var_names = ["ENSX1", "ENSX2", "ENSX3"]
    # ENSX3 carries the backbone target_symbol "Csym" (case-insensitive).
    var_symbols = {"ENSX1": "Asym", "ENSX2": "other", "ENSX3": "CSYM"}
    mapping = map_module_genes(["GENEC"], var_names, var_symbols=var_symbols, backbone_df=bb)
    assert mapping.matched_var_names == ["ENSX3"]
    assert mapping.n_by_symbol == 1


def test_map_module_genes_unmapped():
    bb = _backbone()
    var_names = ["ENSX1", "ENSX2"]
    mapping = map_module_genes(["GENED"], var_names, var_symbols={}, backbone_df=bb)
    assert mapping.matched_var_names == []
    assert mapping.n_by_id == 0
    assert mapping.n_by_symbol == 0


def test_map_module_genes_no_backbone_uses_varname_symbols():
    var_names = ["PRL", "IGFBP1", "FOO"]
    mapping = map_module_genes(["prl", "igfbp1"], var_names, backbone_df=None)
    assert set(mapping.matched_var_names) == {"PRL", "IGFBP1"}
    assert mapping.n_by_symbol == 2


def _make_species(rng: np.random.Generator, n_samples: int, signal: float) -> AnnData:
    """Bulk AnnData: first 4 genes are module genes lifted by `signal`,
    followed by background genes so the within-sample z-score is stable."""
    n_bg = 40
    n_genes = 4 + n_bg
    base = rng.lognormal(mean=3.0, sigma=0.3, size=(n_samples, n_genes))
    base[:, :4] *= np.exp(signal)  # raise the module genes
    symbols = ["GENEA", "GENEB", "GENEC", "GENED"] + [f"BG{i}" for i in range(n_bg)]
    ids = [f"ENSX{i + 1}" for i in range(n_genes)]
    var = pd.DataFrame({"gene_symbol": symbols}, index=ids)
    obs = pd.DataFrame(index=[f"s{i}" for i in range(n_samples)])
    return AnnData(X=base.astype("float32"), obs=obs, var=var)


def test_score_species_pseudobulk_structure_and_floor():
    rng = np.random.default_rng(0)
    adata = _make_species(rng, n_samples=3, signal=1.0)
    bb = pd.DataFrame(
        {
            "source_symbol": ["GENEA", "GENEB", "GENEC", "GENED"],
            "target_gene_id": ["ENSX1", "ENSX2", "ENSX3", "ENSX4"],
            "target_symbol": ["GENEA", "GENEB", "GENEC", "GENED"],
            "tier": [1, 1, 1, 1],
        }
    )
    gene_sets = {"mod_full": ["GENEA", "GENEB", "GENEC", "GENED"], "mod_thin": ["GENEA"]}
    out = score_species_pseudobulk(adata, gene_sets, species="sp", backbone_df=bb)

    assert set(out["score"]) == {"mod_full", "mod_thin"}
    assert (out["species"] == "sp").all()
    full = out[out["score"] == "mod_full"]
    assert len(full) == 3  # one row per sample
    assert (full["n_mapped"] == 4).all()
    assert full["value"].notna().all()

    # mod_thin maps only 1 gene < min_module_genes=3 → NaN scores.
    thin = out[out["score"] == "mod_thin"]
    assert (thin["n_mapped"] == 1).all()
    assert thin["value"].isna().all()


def test_trait_contrast_recovers_planted_signal():
    rng = np.random.default_rng(1)
    frames = []
    # Two trait-positive species with elevated module score, two negative flat.
    for sp, signal in [("pos1", 1.2), ("pos2", 1.0), ("neg1", 0.0), ("neg2", 0.0)]:
        adata = _make_species(rng, n_samples=3, signal=signal)
        bb = pd.DataFrame(
            {
                "source_symbol": ["GENEA", "GENEB", "GENEC", "GENED"],
                "target_gene_id": ["ENSX1", "ENSX2", "ENSX3", "ENSX4"],
                "target_symbol": ["GENEA", "GENEB", "GENEC", "GENED"],
                "tier": [1, 1, 1, 1],
            }
        )
        frames.append(
            score_species_pseudobulk(
                adata, {"mod_full": ["GENEA", "GENEB", "GENEC", "GENED"]}, sp, bb
            )
        )
    scores = pd.concat(frames, ignore_index=True)
    traits = {"pos1": True, "pos2": True, "neg1": False, "neg2": False}

    out = trait_contrast(scores, traits)
    row = out[out["score"] == "mod_full"].iloc[0]
    assert row["n_samples_pos"] == 6
    assert row["n_samples_neg"] == 6
    assert row["delta_pos_minus_neg"] > 0
    assert row["cohens_d"] > 0
    assert row["welch_p"] < 0.05
    assert row["fdr"] <= 1.0


def test_trait_contrast_min_samples_floor():
    scores = pd.DataFrame(
        {
            "species": ["pos1", "neg1"],
            "sample": ["a", "b"],
            "score": ["m", "m"],
            "value": [1.0, 0.0],
        }
    )
    out = trait_contrast(scores, {"pos1": True, "neg1": False}, min_samples_per_arm=2)
    row = out.iloc[0]
    assert np.isnan(row["welch_p"])
    assert np.isnan(row["cohens_d"])
    assert np.isnan(row["fdr"])


def test_trait_contrast_rejects_unknown_species():
    scores = pd.DataFrame({"species": ["x"], "sample": ["a"], "score": ["m"], "value": [1.0]})
    try:
        trait_contrast(scores, {"y": True})
    except ValueError as exc:
        assert "no entry for species" in str(exc)
    else:  # pragma: no cover
        raise AssertionError("expected ValueError")
