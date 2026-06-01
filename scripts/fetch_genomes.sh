#!/usr/bin/env bash
# Fetch the genome FASTAs for the Q4.5 cross-species trigger-element test.
#
# Downloads, decompresses and indexes (.fai) the four genomes of the minimal
# convergence set. Re-running is cheap: a genome already present (indexed) is
# skipped. Indexing uses pysam (pip-native) so no samtools binary is required.
#
#   Human            hg19            UCSC      -> results/raw/reference/hg19.fa
#   Mouse            mm10            UCSC      -> results/raw/genomes/mouse.fa
#   Ground squirrel  speTri2         UCSC      -> results/raw/genomes/ground_squirrel.fa
#   Carollia bat     GCA_056371365.1 NCBI/Bat1K-> results/raw/genomes/bat_carollia.fa
#
# The bat is the ACTUAL Carollia perspicillata assembly (mCarPer1.2, the
# chromosome-level Bat1K build) — NOT the Myotis lucifugus reference the RNA
# was aligned to. Carollia is the menstruating/spontaneous species under test.
set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REF_DIR="$ROOT/results/raw/reference"
GEN_DIR="$ROOT/results/raw/genomes"
mkdir -p "$REF_DIR" "$GEN_DIR"

UCSC="https://hgdownload.cse.ucsc.edu/goldenPath"
NCBI="https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/056/371/365/GCA_056371365.1_mCarPer1.2"

# label  dest_fasta  url
GENOMES=(
  "hg19|$REF_DIR/hg19.fa|$UCSC/hg19/bigZips/hg19.fa.gz"
  "mm10|$GEN_DIR/mouse.fa|$UCSC/mm10/bigZips/mm10.fa.gz"
  "speTri2|$GEN_DIR/ground_squirrel.fa|$UCSC/speTri2/bigZips/speTri2.fa.gz"
  "mCarPer1.2|$GEN_DIR/bat_carollia.fa|$NCBI/GCA_056371365.1_mCarPer1.2_genomic.fna.gz"
)

faidx() { python -c "import pysam,sys; pysam.faidx(sys.argv[1])" "$1"; }

for entry in "${GENOMES[@]}"; do
  IFS='|' read -r label fa url <<<"$entry"
  if [[ -f "$fa.fai" ]]; then
    echo "[skip] $label already present ($fa)"
    continue
  fi
  echo "[get ] $label <- $url"
  curl -fL --retry 3 -o "$fa.gz" "$url"
  echo "[gunzip] $fa.gz"
  gunzip -f "$fa.gz"
  echo "[faidx] $fa"
  faidx "$fa"
  echo "[done] $label -> $fa ($(du -h "$fa" | cut -f1))"
done

echo "All genomes ready:"
ls -lh "$REF_DIR/hg19.fa" "$GEN_DIR"/*.fa
