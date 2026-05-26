"""Smoke tests against real downloaded data.

Skipped by default in CI. Run locally after the pipeline has produced
real artifacts in ``results/``::

    pytest -m real_data
"""

from __future__ import annotations

from pathlib import Path

import pytest

REAL_DATA = pytest.mark.real_data
PROJECT_ROOT = Path(__file__).resolve().parents[1]
RESULTS = PROJECT_ROOT / "results"


@REAL_DATA
def test_integrated_stromal_exists():
    p = RESULTS / "integrated" / "stromal_harmony.h5ad"
    if not p.exists():
        pytest.skip(f"{p} missing; run `wombat integrate --mode stromal` first")
    import anndata as ad

    a = ad.read_h5ad(p, backed="r")
    try:
        assert a.n_obs > 0
        assert "species" in a.obs.columns
    finally:
        a.file.close()


@REAL_DATA
def test_integrated_includes_mouse():
    """Q1 exit criterion: mouse stromal cells must be in the joint object."""
    p = RESULTS / "integrated" / "stromal_harmony.h5ad"
    if not p.exists():
        pytest.skip("integrated h5ad missing")
    import anndata as ad

    a = ad.read_h5ad(p, backed="r")
    try:
        species = set(a.obs["species"].astype(str).unique())
    finally:
        a.file.close()
    assert "mouse" in species, f"Mouse not yet integrated. species={species}"


@REAL_DATA
def test_mouse_stromal_recall():
    """Q2.1/Q2.4 regression guard: ≥60% of UE_DSC mouse cells must
    survive annotation + stromal subset. GSE226417 is "Uterine Epithelial
    AND Decidual Stromal" by design — the ~33% non-stromal fraction is
    largely real epithelial signal (Muc1+Krt18+Epcam), not annotation
    failure. The Q2.1 80% aspiration was based on the incorrect
    assumption that the dataset was pure stromal; replaced in Q2.4 with
    the hierarchical lineage gate (architectural improvement) plus this
    honest 60% floor.

    Baseline before species_overrides: 11,484 / 23,471 = 48.9%
    With species_overrides (Q2.1):     15,662 / 23,471 = 66.7%
    With hierarchical lineage (Q2.4):  15,662 / 23,471 = 66.7%
    """
    p = RESULTS / "integrated" / "stromal_harmony.h5ad"
    if not p.exists():
        pytest.skip("integrated h5ad missing")
    import anndata as ad

    qc_path = RESULTS / "qc" / "GSE226417.h5ad"
    if not qc_path.exists():
        pytest.skip("GSE226417 QC h5ad missing")

    a = ad.read_h5ad(p, backed="r")
    try:
        if "dataset" not in a.obs.columns:
            pytest.skip("dataset column not in integrated obs")
        mouse_stromal = int((a.obs["dataset"].astype(str) == "GSE226417").sum())
    finally:
        a.file.close()

    qc = ad.read_h5ad(qc_path, backed="r")
    try:
        ue_dsc_total = qc.n_obs
    finally:
        qc.file.close()

    recall = mouse_stromal / ue_dsc_total
    assert recall >= 0.60, (
        f"Mouse stromal recall regressed: {mouse_stromal}/{ue_dsc_total} "
        f"= {recall:.1%} (Q2.1 floor 60%; aspirational 80%)"
    )


@REAL_DATA
def test_manifest_has_checksums():
    p = RESULTS / "reports" / "manifest.csv"
    if not p.exists():
        pytest.skip("manifest.csv missing; run `wombat generate-reports`")
    import pandas as pd

    df = pd.read_csv(p)
    assert "sha256" in df.columns
    assert df["sha256"].notna().all(), "manifest has rows with null sha256"
    assert (df["sha256"].str.len() == 64).all(), "non-sha256 hash present"


@REAL_DATA
def test_coverage_matches_disk():
    """Coverage report's `integrated` column must reflect actual h5ad contents."""
    import anndata as ad
    import pandas as pd

    coverage = RESULTS / "reports" / "coverage.md"
    integrated = RESULTS / "integrated" / "stromal_harmony.h5ad"
    if not coverage.exists() or not integrated.exists():
        pytest.skip("coverage.md or integrated h5ad missing")

    a = ad.read_h5ad(integrated, backed="r")
    try:
        actual = (
            set(a.obs["dataset"].astype(str).unique()) if "dataset" in a.obs.columns else set()
        )
    finally:
        a.file.close()

    text = coverage.read_text()
    # Every claimed-integrated accession in the report must actually appear in the h5ad.
    for line in text.splitlines():
        if "True" in line and "scRNA" in line:
            acc = line.split("|")[1].strip() if "|" in line else ""
            if acc.startswith(("GSE", "E-MTAB")):
                # only enforce if the report shows integrated=True for it
                cells = [c.strip() for c in line.split("|")]
                if len(cells) >= 7 and cells[5] == "True":
                    assert acc in actual, (
                        f"coverage.md marks {acc} as integrated but it's not in "
                        f"the integrated h5ad ({sorted(actual)})"
                    )
    _ = pd  # keep import for future expansion


@REAL_DATA
def test_protected_core_markers_survive():
    """Pre-Q3 gate item A/B regression: the protected core decidual
    panel must survive HVG selection and remain in the integrated
    h5ad's ``var_names``.

    Source of truth: ``configs/markers.yaml::protected_core``. Any
    gene listed there that is **absent** from the joint upstream
    var set (so the HVG carveout cannot recover it) is reported as
    a skip with the explanation, not a failure — that case requires
    an orthology or upstream-QC fix, not an integration-code fix.
    """
    p = RESULTS / "integrated" / "stromal_harmony.h5ad"
    if not p.exists():
        pytest.skip("integrated h5ad missing")

    from wombat.config import load_config

    markers = load_config("markers")
    core = markers.get("protected_core") if isinstance(markers, dict) else None
    if not core:
        pytest.skip("protected_core not defined in configs/markers.yaml")

    import anndata as ad

    a = ad.read_h5ad(p, backed="r")
    try:
        var_set = {str(v) for v in a.var_names}
    finally:
        a.file.close()

    # Recoverable = present in at least one upstream processed h5ad
    # (so the HVG carveout in integrate.py could have force-included it).
    processed_dir = RESULTS / "processed"
    upstream: set[str] = set()
    if processed_dir.exists():
        for ph in processed_dir.glob("*.h5ad"):
            ph_a = ad.read_h5ad(ph, backed="r")
            try:
                upstream |= {str(v).upper() for v in ph_a.var_names}
            finally:
                ph_a.file.close()

    recoverable = [g for g in core if g.upper() in upstream] if upstream else core
    not_recoverable = sorted(set(core) - set(recoverable))
    missing_recoverable = sorted(g for g in recoverable if g not in var_set)

    if not_recoverable and not missing_recoverable:
        pytest.skip(
            "All recoverable protected-core markers present; the following "
            f"are upstream-absent (not an integration bug): {not_recoverable}"
        )

    assert not missing_recoverable, (
        "Protected-core markers dropped from integrated var set despite "
        f"being present upstream — HVG carveout is broken. Missing: "
        f"{missing_recoverable}. Upstream-absent (out of scope): "
        f"{not_recoverable}."
    )
