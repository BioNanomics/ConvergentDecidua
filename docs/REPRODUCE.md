# Reproducing the ConvergentDecidua Human–Mouse Atlas

This guide walks an external reviewer from a clean machine to a reproduced
human–mouse stromal decidualization atlas with checksums.

> **Status:** Q1 of the 12-month plan — mouse scRNA ingest is the critical
> path. If any step here fails, please open an issue with the failing command
> and full output. That is the bar we hold ourselves to.

---

## 1. Prerequisites

- **Python** 3.9–3.12 (project tested on 3.9 and 3.11).
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

```bash
pytest -q                    # unit tests (no real data needed)
pytest -m real_data -q       # smoke tests against results/ (skipped without data)
```

Key real-data assertions:

- `results/integrated/stromal_harmony.h5ad` exists and contains both
  `human` and `mouse` in `.obs.species`.
- `results/reports/manifest.csv` has a 64-char sha256 for every row.
- `results/reports/coverage.md` only marks an accession as `integrated=True`
  when its cells actually appear in the integrated h5ad.

---

## 5. Known limitations (Q1)

- **scATAC (GSE183771)** is in scope for Q2; ingest scaffolding exists but
  joint analysis is not yet runnable end-to-end.
- **E-MTAB-11491** (mouse, 644 per-sample files) is a Q2 stretch fallback,
  not a Q1 dependency. GSE226417 alone supplies the mouse stromal cells.
- **PRL** has no 1:1 mouse ortholog; the mouse decidual-prolactin family
  (`Prl8a2`, `Prl3c1`, `Prl3d1`) will be added via per-species score
  overrides in Q3.
