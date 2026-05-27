"""Tests for the species-driven Ensembl ortholog dispatch.

These tests check that `src.orthologs.ensembl` looks up dataset names,
attribute prefixes, and Compara species identifiers from
`configs/species.yaml` rather than from any hard-coded table. No
network calls are made — the BioMart HTTP layer is mocked.
"""

from __future__ import annotations

from unittest.mock import patch

import pytest

from src.orthologs.ensembl import (
    _ortholog_attrs,
    _parse_biomart_response,
    _species_field,
    _species_index,
    fetch_ensembl_orthologs,
)


def test_species_index_has_required_fields():
    idx = _species_index()
    # Sanity: at minimum the Tier A pair must be there with all three
    # generalization-relevant fields.
    for name in ("human", "mouse"):
        assert name in idx, f"{name} missing from species.yaml"
        for field in ("ensembl_dataset", "ensembl_prefix", "ensembl_species"):
            assert idx[name].get(field), f"{name}.{field} missing in species.yaml"


def test_species_field_resolution_for_macaque():
    """Macaque was switched from rhesus → cynomolgus 2026-05-27. Lookups
    must reflect the cynomolgus identifiers, not the old rhesus ones."""
    assert _species_field("macaque", "ensembl_dataset") == "mfascicularis_gene_ensembl"
    assert _species_field("macaque", "ensembl_prefix") == "mfascicularis"
    assert _species_field("macaque", "ensembl_species") == "macaca_fascicularis"


def test_species_field_unknown_species_raises():
    with pytest.raises(ValueError, match="Unknown species"):
        _species_field("nonexistent_species", "ensembl_dataset")


def test_ortholog_attrs_uses_target_prefix():
    attrs = _ortholog_attrs("mfascicularis")
    assert attrs["homolog_ensembl_gene"] == "mfascicularis_homolog_ensembl_gene"
    assert attrs["homolog_orthology_type"] == "mfascicularis_homolog_orthology_type"
    assert attrs["homolog_orthology_confidence"] == ("mfascicularis_homolog_orthology_confidence")


def test_parse_biomart_response_finds_target_columns_for_any_species():
    """The parser must locate target stable-ID/name columns without
    knowing the target's display label (mouse, macaque, baboon, ...)."""
    tsv = (
        "Gene stable ID\tGene name\tBaboon gene stable ID\tBaboon gene name\t"
        "Baboon homology type\tBaboon orthology confidence [0 low, 1 high]\n"
        "ENSG00000001\tHAND2\tENSPANG000099\thand2\tortholog_one2one\t1\n"
        "ENSG00000002\tFOXO1\tENSPANG000088\tfoxo1\tortholog_one2one\t1\n"
    )
    tbl = _parse_biomart_response(tsv, "human", "baboon")
    assert tbl.num_rows == 2
    assert set(tbl.column_names) == {
        "human_gene_id",
        "human_symbol",
        "baboon_gene_id",
        "baboon_symbol",
        "orthology_type",
        "confidence",
    }
    assert tbl.column("baboon_gene_id").to_pylist() == ["ENSPANG000099", "ENSPANG000088"]


def test_fetch_uses_species_yaml_dataset(monkeypatch):
    """fetch_ensembl_orthologs must build the XML query using the
    dataset name from species.yaml, not from a hard-coded dict."""
    captured: dict = {}

    class _FakeResp:
        text = (
            "Gene stable ID\tGene name\tMouse gene stable ID\tMouse gene name\t"
            "Mouse homology type\tMouse orthology confidence [0 low, 1 high]\n"
            "ENSG00000001\tHAND2\tENSMUSG000099\thand2\tortholog_one2one\t1\n"
        )
        headers = {"content-type": "text/tab-separated-values"}

        def raise_for_status(self):
            return None

    def _fake_get(url, params=None, timeout=None):
        captured["url"] = url
        captured["xml"] = params["query"] if params else ""
        return _FakeResp()

    with patch("requests.get", _fake_get):
        tbl = fetch_ensembl_orthologs("human", "mouse", cache_dir=None)

    assert tbl.num_rows == 1
    # Both lookups should appear in the XML — proving the values came
    # from species.yaml, not a hard-coded fallback.
    assert "hsapiens_gene_ensembl" in captured["xml"]
    assert "mmusculus_homolog_ensembl_gene" in captured["xml"]
