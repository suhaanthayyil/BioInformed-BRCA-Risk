#!/bin/bash
# Download raw data for the breast cancer prediction project.
# Run this script from the data/raw/ directory.
#
# Covers all four cohorts:
#   - TCGA-BRCA (UCSC Xena HiSeqV2, clinical matrix, TCGA-CDR survival)
#   - GSE96058 / SCAN-B (GEO expression + series matrices)
#   - METABRIC (cBioPortal datahub)
#   - GSE20685 (GEO series matrix + GPL570 annotation)
# Gene sets: see download_gene_sets.sh. Then run scripts/ingest_raw.py.
#
# SHA256 values are cross-checked against the *.metadata.json files in
# data/processed/. The script fails loudly on any mismatch.

set -euo pipefail

verify() { echo "$1  $2" | shasum -a 256 -c -; }

echo "=== [0/3] TCGA-BRCA (UCSC Xena, ~65MB) ==="
mkdir -p tcga gse96058
XENA="https://tcga.xenahubs.net/download"
[ -f tcga/HiSeqV2.gz ] || curl -fL -o tcga/HiSeqV2.gz "$XENA/TCGA.BRCA.sampleMap/HiSeqV2.gz"
[ -f tcga/BRCA_clinicalMatrix ] || curl -fL -o tcga/BRCA_clinicalMatrix "$XENA/TCGA.BRCA.sampleMap/BRCA_clinicalMatrix"
[ -f tcga/BRCA_survival.txt ] || curl -fL -o tcga/BRCA_survival.txt "$XENA/survival/BRCA_survival.txt"
verify "263bf67245cc4b9062676583c0ff0306f08471a26aafd5504037e1da22133746" tcga/HiSeqV2.gz
verify "39eb3be0fb86e6a577bd2cc01502a7fa5a271e1e1cba294e9dc644ad99580d7f" tcga/BRCA_clinicalMatrix
verify "6c9e18ae9c2b8bacc48cd5a60cde97e8880ae2e84e48df437a0378b17e659b73" tcga/BRCA_survival.txt

echo "=== [1/3] GSE96058 / SCAN-B gene expression (~565MB compressed) + clinical ==="
GEO96="https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/matrix"
for f in GSE96058-GPL11154_series_matrix.txt.gz GSE96058-GPL18573_series_matrix.txt.gz; do
  [ -f "gse96058/$f" ] || curl -fL -o "gse96058/$f" "$GEO96/$f"
done
verify "47e4e488911f688a9c9ac7a0bb1223a22259ca1dcae2a1a7501ae5b5e55794ce" gse96058/GSE96058-GPL11154_series_matrix.txt.gz
verify "d7a77604b78dc01235ebee55e39aceb4ff24ecae4fd540e3fa2a773ddb8be091" gse96058/GSE96058-GPL18573_series_matrix.txt.gz
GSE96058_GZ="GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz"
if [ ! -f "$GSE96058_GZ" ] && [ ! -f "GSE96058_gene_expression.csv" ]; then
  curl -fL -o "$GSE96058_GZ" \
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/${GSE96058_GZ}"
fi
[ -f "$GSE96058_GZ" ] && verify "3c717baf7960e1f1477a72744a399a42fabba113c9e7ecb3d55997574d7e9732" "$GSE96058_GZ"

echo "=== [2/3] METABRIC (cBioPortal datahub, ~150MB) ==="
METABRIC_TAR="brca_metabric.tar.gz"
if [ ! -f "$METABRIC_TAR" ] && [ ! -d metabric ]; then
  curl -fL -o "$METABRIC_TAR" "https://datahub.assets.cbioportal.org/brca_metabric.tar.gz"
fi
if [ -f "$METABRIC_TAR" ]; then
  verify "6d4683477d6b37a2d7edbedc0df610f67bc456f99e5e1bef6219f37b633a55f7" "$METABRIC_TAR"
  mkdir -p metabric && tar -xzf "$METABRIC_TAR" -C metabric --strip-components=1
fi

echo "=== [3/3] GSE20685 series matrix + GPL570 annotation (~145MB) ==="
GSE20685_GZ="GSE20685_series_matrix.txt.gz"
GPL570_GZ="GPL570.annot.gz"
[ -f "$GSE20685_GZ" ] || curl -fL -o "$GSE20685_GZ" \
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE20nnn/GSE20685/matrix/${GSE20685_GZ}"
[ -f "$GPL570_GZ" ] || curl -fL -o "$GPL570_GZ" \
  "https://ftp.ncbi.nlm.nih.gov/geo/platforms/GPL570/annot/${GPL570_GZ}"
verify "e818a4d5834e20bbcf01de515caf856c394abb36a02b9cccc95be05c22d0a279" "$GSE20685_GZ"
verify "d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394" "$GPL570_GZ"

echo "=== Done. Next: bash download_gene_sets.sh, then see data/processed/REGENERATION.md ==="
