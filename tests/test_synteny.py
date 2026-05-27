"""Unit tests for src/orthologs/synteny.py.

Mocks the Ensembl REST endpoint; no network calls.
"""

from __future__ import annotations

from pathlib import Path

import pyarrow.parquet as pq

from src.orthologs.synteny import (
    _parse_homology,
    check_synteny,
    run_synteny_check,
)

# Canonical IGFBP1 human→mouse payload shape (trimmed to the fields we read).
_FAKE_IGFBP1_PAYLOAD = {
    "data": [
        {
            "id": "ENSG00000146678",
            "homologies": [
                {
                    "type": "ortholog_one2one",
                    "method_link_type": "ENSEMBL_ORTHOLOGUES",
                    "dn_ds": 0.221,
                    "source": {
                        "id": "ENSG00000146678",
                        "species": "homo_sapiens",
                        "perc_id": 85.4,
                    },
                    "target": {
                        "id": "ENSMUSG00000020429",
                        "species": "mus_musculus",
                        "symbol": "Igfbp1",
                        "perc_id": 84.1,
                    },
                }
            ],
        }
    ]
}


def test_parse_one2one():
    rows = _parse_homology(_FAKE_IGFBP1_PAYLOAD, "human", "IGFBP1", "mouse")
    assert len(rows) == 1
    r = rows[0]
    assert r["anchor_symbol"] == "IGFBP1"
    assert r["target_symbol"] == "Igfbp1"
    assert r["orthology_type"] == "ortholog_one2one"
    assert bool(r["alignment_present"]) is True
    assert r["perc_id_source"] == 85.4
    assert r["dn_ds"] == 0.221


def test_parse_empty_payload():
    rows = _parse_homology(None, "human", "MISSING", "mouse")
    assert len(rows) == 1
    assert rows[0]["orthology_type"] == "no_data"
    assert bool(rows[0]["alignment_present"]) is False


def test_check_synteny_batches(monkeypatch):
    def fake_fetch(anchor, symbol, target, **kwargs):
        if symbol == "IGFBP1":
            return _FAKE_IGFBP1_PAYLOAD
        return None  # simulate "unknown symbol"

    monkeypatch.setattr("src.orthologs.synteny._fetch_homology", fake_fetch)
    monkeypatch.setattr("src.orthologs.synteny.time.sleep", lambda _s: None)

    table = check_synteny(["IGFBP1", "BOGUS"], ["mouse"])
    df = table.to_pandas()
    assert len(df) == 2
    assert df.loc[df["anchor_symbol"] == "IGFBP1", "alignment_present"].iat[0]
    assert not df.loc[df["anchor_symbol"] == "BOGUS", "alignment_present"].iat[0]


def test_run_writes_parquet(tmp_path: Path, monkeypatch):
    monkeypatch.setattr(
        "src.orthologs.synteny._fetch_homology",
        lambda *a, **kw: _FAKE_IGFBP1_PAYLOAD,
    )
    monkeypatch.setattr("src.orthologs.synteny.time.sleep", lambda _s: None)

    out = tmp_path / "synteny.parquet"
    table = run_synteny_check(out, symbols=["IGFBP1"], target_species=["mouse"])
    assert out.exists()
    roundtrip = pq.read_table(out).to_pandas()
    assert len(roundtrip) == 1
    assert roundtrip.iloc[0]["target_symbol"] == "Igfbp1"
    assert table.num_rows == 1


def test_rest_error_recorded(monkeypatch):
    def boom(*a, **kw):
        raise RuntimeError("connection refused")

    monkeypatch.setattr("src.orthologs.synteny._fetch_homology", boom)
    monkeypatch.setattr("src.orthologs.synteny.time.sleep", lambda _s: None)

    table = check_synteny(["IGFBP1"], ["mouse"])
    df = table.to_pandas()
    assert len(df) == 1
    assert df.iloc[0]["orthology_type"].startswith("rest_error:")
    assert bool(df.iloc[0]["alignment_present"]) is False


def test_cli_synteny_check_help():
    """Smoke-test the wired CLI command."""
    from click.testing import CliRunner

    from wombat.cli import cli

    runner = CliRunner()
    result = runner.invoke(cli, ["orthologs", "synteny-check", "--help"])
    assert result.exit_code == 0
    assert "Ensembl REST" in result.output
