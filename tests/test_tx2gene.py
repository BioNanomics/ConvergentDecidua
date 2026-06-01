"""Unit tests for src/ingest/tx2gene.py.

All BioMart/REST access is monkeypatched; no network calls.
"""

from __future__ import annotations

import anndata as ad
import numpy as np
import pyarrow as pa
import pytest
import scipy.sparse as sp

from src.ingest import tx2gene as t2g


def test_build_xml_without_chromosome():
    xml = t2g._build_tx2gene_xml("mmusculus_gene_ensembl")
    assert 'name="mmusculus_gene_ensembl"' in xml
    assert "chromosome_name" not in xml
    assert 'name="ensembl_transcript_id"' in xml


def test_build_xml_with_chromosome_filter():
    xml = t2g._build_tx2gene_xml("mmusculus_gene_ensembl", chromosomes="3")
    assert '<Filter name="chromosome_name" value="3" />' in xml


def test_group_chromosome_chunks_splits_main_and_batches_scaffolds():
    regions = ["1", "2", "X", "Y", "MT"] + [f"SCAF{i:04d}AB" for i in range(120)]
    chunks = t2g._group_chromosome_chunks(regions)
    # 5 karyotype regions emitted one-per-chunk, 120 scaffolds -> ceil(120/50)=3
    assert chunks[:5] == ["1", "2", "X", "Y", "MT"]
    scaffold_chunks = chunks[5:]
    assert len(scaffold_chunks) == 3
    assert scaffold_chunks[0].count(",") + 1 == t2g._SCAFFOLD_BATCH
    # Every scaffold appears exactly once across the batched chunks.
    flattened = ",".join(scaffold_chunks).split(",")
    assert len(flattened) == 120
    assert len(set(flattened)) == 120


def test_parse_tsv_drops_rows_missing_ids():
    text = (
        "Transcript stable ID\tGene stable ID\tGene name\n"
        "ENSMUST1\tENSMUSG1\tFoo\n"
        "\tENSMUSG2\tBar\n"  # missing transcript -> dropped
        "ENSMUST3\t\tBaz\n"  # missing gene -> dropped
        "ENSMUST4\tENSMUSG4\t\n"  # missing symbol -> kept
    )
    tbl = t2g._parse_tx2gene_tsv(text)
    df = tbl.to_pandas()
    assert list(df["transcript_id"]) == ["ENSMUST1", "ENSMUST4"]
    assert list(df["gene_id"]) == ["ENSMUSG1", "ENSMUSG4"]
    assert list(df["gene_symbol"]) == ["Foo", ""]


def test_parse_tsv_rejects_unexpected_columns():
    with pytest.raises(RuntimeError, match="Unexpected BioMart columns"):
        t2g._parse_tx2gene_tsv("colA\tcolB\nx\ty\n")


def test_strip_version():
    assert t2g._strip_version(["ENSMUST1.3", "ENSMUST2", "X.10"]) == [
        "ENSMUST1",
        "ENSMUST2",
        "X",
    ]


def test_fetch_chunk_with_retry_recovers(monkeypatch):
    calls = {"n": 0}

    def flaky_then_ok(_xml, _label):
        calls["n"] += 1
        if calls["n"] < 3:
            raise RuntimeError("all mirrors failed")
        return "ok-body", []

    monkeypatch.setattr(t2g, "_fetch_biomart_tsv", flaky_then_ok)
    monkeypatch.setattr("time.sleep", lambda _s: None)

    text, flaky = t2g._fetch_chunk_with_retry("<xml/>", "mouse:3")
    assert text == "ok-body"
    assert calls["n"] == 3
    # The two transient failures are recorded for transparency.
    assert len(flaky) == 2


def test_fetch_chunk_with_retry_gives_up(monkeypatch):
    def always_fail(_xml, _label):
        raise RuntimeError("all mirrors failed")

    monkeypatch.setattr(t2g, "_fetch_biomart_tsv", always_fail)
    monkeypatch.setattr("time.sleep", lambda _s: None)

    with pytest.raises(RuntimeError, match="all mirrors failed"):
        t2g._fetch_chunk_with_retry("<xml/>", "mouse:3")


def test_fetch_tx2gene_chunked_concatenates_and_dedups(monkeypatch):
    monkeypatch.setattr(t2g, "_fetch_chromosome_names", lambda _sp: ["1", "2"])

    bodies = {
        "1": (
            "Transcript stable ID\tGene stable ID\tGene name\n"
            "ENSMUST1\tENSMUSG1\tFoo\n"
            "ENSMUST2\tENSMUSG1\tFoo\n"
        ),
        "2": (
            "Transcript stable ID\tGene stable ID\tGene name\n"
            "ENSMUST3\tENSMUSG2\tBar\n"
            "ENSMUST1\tENSMUSG1\tFoo\n"  # duplicate transcript across chunks
        ),
    }

    def fake_fetch(_xml, label):
        chrom = label.split(":", 1)[1]
        return bodies[chrom], []

    monkeypatch.setattr(t2g, "_fetch_biomart_tsv", fake_fetch)
    monkeypatch.setattr("time.sleep", lambda _s: None)

    table, flaky = t2g._fetch_tx2gene_chunked("mmusculus_gene_ensembl", "mus_musculus", "mouse")
    df = table.to_pandas()
    assert flaky == []
    # 3 unique transcripts after de-dup (ENSMUST1 appeared twice).
    assert sorted(df["transcript_id"]) == ["ENSMUST1", "ENSMUST2", "ENSMUST3"]
    assert df["transcript_id"].nunique() == len(df)


def _tx2gene_table() -> pa.Table:
    return pa.table(
        {
            "transcript_id": ["ENSMUST1", "ENSMUST2", "ENSMUST3"],
            "gene_id": ["ENSMUSG1", "ENSMUSG1", "ENSMUSG2"],
            "gene_symbol": ["Foo", "Foo", "Bar"],
        }
    )


def test_aggregate_to_genes_sums_transcripts():
    # Two transcripts of gene ENSMUSG1 + one of ENSMUSG2; include an
    # unmapped transcript that must be dropped.
    X = np.array([[1.0, 2.0, 5.0, 9.0], [0.0, 4.0, 6.0, 1.0]])
    adata = ad.AnnData(
        X=X,
        var={"_": [0, 0, 0, 0]},
    )
    adata.var_names = ["ENSMUST1.2", "ENSMUST2", "ENSMUST3", "ENSMUST_UNKNOWN"]

    out = t2g.aggregate_to_genes(adata, _tx2gene_table())
    assert list(out.var_names) == ["ENSMUSG1", "ENSMUSG2"]
    # ENSMUSG1 = ENSMUST1 + ENSMUST2 per row; ENSMUSG2 = ENSMUST3.
    np.testing.assert_array_equal(np.asarray(out.X), np.array([[3.0, 5.0], [4.0, 6.0]]))
    assert list(out.var["gene_symbol"]) == ["Foo", "Bar"]


def test_aggregate_to_genes_sparse_matches_dense():
    X = sp.csr_matrix(np.array([[1.0, 2.0, 5.0], [0.0, 4.0, 6.0]]))
    adata = ad.AnnData(X=X)
    adata.var_names = ["ENSMUST1", "ENSMUST2", "ENSMUST3"]

    out = t2g.aggregate_to_genes(adata, _tx2gene_table())
    dense = np.asarray(out.X.todense()) if sp.issparse(out.X) else np.asarray(out.X)
    np.testing.assert_array_equal(dense, np.array([[3.0, 5.0], [4.0, 6.0]]))


def test_aggregate_to_genes_raises_when_nothing_maps():
    adata = ad.AnnData(X=np.array([[1.0, 2.0]]))
    adata.var_names = ["NOPE1", "NOPE2"]
    with pytest.raises(RuntimeError, match="No transcripts.*mapped"):
        t2g.aggregate_to_genes(adata, _tx2gene_table())
