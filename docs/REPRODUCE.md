# Reproducing the ConvergentDecidua Human–Mouse Atlas

This guide walks an external reviewer from a clean machine to a reproduced
human–mouse stromal decidualization atlas with checksums.

> **Status:** Q2 closed (May 2026). Q3 is blocked on a pre-Q3
> acceptance gate — see [`../PLAN.md`](../PLAN.md) and
> [`q2_closeout.md`](q2_closeout.md) for what the current pipeline can
> and cannot support scientifically. If any step here fails, please
> open an issue with the failing command and full output.

---

## 1. Prerequisites

- **Python** 3.11 or 3.12 (`requires-python = ">=3.11,<3.13"`).
- **R** ≥ 4.1 with `Seurat` and `Matrix` packages — only needed for the
  GSE226417 mouse ingest. Skip if you do not need to re-ingest the mouse
  RData. The Docker image at `docker/Dockerfile` provides this pre-built.
- **~30 GB free disk** for raw + processed artifacts.
- **Network access** to GEO and Ensembl Compara FTP.

### Option A — local install

```bash
git clone https://github.com/BioNanomics/ConvergentDecidua
cd ConvergentDecidua
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev,ingest,workflow]"
wombat --help
```

### Option B — Docker (includes R + Seurat)

```bash
# Build from the repo root so the Dockerfile can COPY src/, wombat/, etc.
docker build -t convergent-decidua -f docker/Dockerfile .
docker run --rm -v "$PWD/results:/app/results" convergent-decidua --help
```

---

## 2. Pipeline (CLI form)

Every step has a Snakemake equivalent under `workflows/`; the CLI form below
is the source of truth.

```bash
# 1. Configs + ortholog backbone (no large downloads beyond Ensembl Compara)
wombat validate-config
wombat build-registry
wombat orthologs build

# 2. Fetch + convert each dataset
wombat fetch --dataset GSE111976           # human scRNA
wombat fetch --dataset GSE127918           # human scRNA
wombat fetch --dataset GSE226429           # mouse bulk
wombat fetch --dataset GSE226417           # mouse scRNA (REQUIRES R + Seurat)

# 3. QC
wombat qc --species human
wombat qc --species mouse

# 4. Cross-species stromal integration
wombat integrate --mode stromal --method harmony

# 5. Decidualization scoring
wombat score-decidua

# 6. Reports + manifest with sha256 checksums
wombat generate-reports
```

After step 6, `results/reports/manifest.csv` contains a `sha256` for every
artifact in `results/`. Compare against the published manifest in the release
to confirm bit-identical reproduction.

---

## 3. Mouse scRNA (GSE226417) — the R bridge

GSE226417 ships as four Seurat `.RData.gz` dumps totaling ~10 GB. The dataset
registry (`configs/datasets.yaml`) declares:

```yaml
ingest:
  format: seurat_rdata
  include: [UE_DSC]
```

This restricts conversion to the decidual-stromal subset (~1.5 GB), which is
the population we need for the cross-species stromal atlas. To include other
tissues, edit `include` (e.g. `[UE_DSC, UE_EC]`) or remove it to convert all
four files.

The bridge:

1. `wombat fetch --dataset GSE226417` downloads the `.RData.gz` files into
   `results/raw/GSE226417/`.
2. `src/ingest/anndata_writer.py` detects the `seurat_rdata` ingest format and
   shells out via `Rscript scripts/rdata_to_h5ad.R`.
3. The R script emits a 10X-style MTX directory plus `obs.csv` under
   `results/raw/GSE226417/_rdata_export/<file_stem>/`.
4. Python loads those files into a standardized AnnData and writes
   `results/processed/GSE226417.h5ad`.

If `Rscript` is not on PATH the Python layer raises a clear error pointing
back to this document.

---

## 4. Verifying a reproduction

### What CI verifies on every push

GitHub Actions (`.github/workflows/ci.yml`) runs four jobs on every
push and pull request. These are **code-quality and workflow-syntax**
checks; they do **not** download data or run the biological pipeline.

| Job | What it checks | Why it is not "biological reproduction" |
|---|---|---|
| `lint` | `ruff check .` + `ruff format --check .` | Code style only |
| `test` | `pytest -q` (default markers) | Unit tests on synthetic fixtures; `real_data` marker is **deselected** in `pyproject.toml` |
| `validate-configs` | `wombat validate-config` | YAML schema only |
| `validate-workflow` | `snakemake -n --snakefile workflows/Snakefile --forceall` | DAG resolution only — no rules execute |

A green CI badge means "the code compiles, the unit tests pass on
synthetic data, the configs parse, and the workflow DAG resolves."
**It does not mean the biological pipeline has been re-run.**

### What you (the reviewer) must run locally for real-data reproduction

After running the full pipeline in §2:

```bash
pytest -m real_data -q          # smoke tests against results/
                                # (skipped in CI; required here)
```

Key real-data assertions:

- `results/integrated/stromal_cross_species.h5ad` exists and contains
  both `human` and `mouse` in `.obs.species`
  (`stromal_harmony.h5ad` is a symlink for back-compat).
- Mouse stromal recall ≥ 60 % on GSE226417 UE_DSC
  (`test_mouse_stromal_recall`).
- `results/reports/manifest.csv` has a 64-char sha256 for every row.
- `results/reports/coverage.md` only marks an accession as
  `integrated=True` when its cells actually appear in the integrated
  h5ad.

If any of those fail, please open an issue with the command and the
full output.

---

## 5. Known limitations (post-Q2)

- **Cross-species mixing (LISI ≈ 1.00).** Harmony does not produce a
  shared embedding on the current joint feature space. Pre-Q3 gate
  items A + E address this.
- **5 / 8 canonical decidual markers dropped by HVG selection**
  (PGR, HAND2, WNT4, PRL, LEFTY2). Pre-Q3 gate item A introduces a
  `protected_core` carveout. Do not interpret module-level conserved/
  divergent claims from the current integrated h5ad.
- **Orthology backbone is not externally validated** (g:Profiler
  confirmed 0 / 16 168 Tier 1 mappings). Pre-Q3 gate item C requires
  a per-gene spot-check memo for the protected core before any
  comparative-biology claim.
- **scATAC (GSE183771)** preprocessing primitives exist
  (`src/qc/scatac.py`) but fetch + processing are deferred to Q3
  stretch.
- **PRL has no 1:1 mouse ortholog.** `species_overrides` handles
  scoring; PRL/LEFTY2 are explicitly exploratory and will not carry
  year-one claims (see `marker_recovery_plan.md`).
- **Python 3.11 is the floor** (`requires-python = ">=3.11,<3.13"`).
  3.9/3.10 are no longer supported.
