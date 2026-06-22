#!/bin/bash
# Download raw data for the breast cancer prediction project.
# Run this script from the data/raw/ directory.
#
# Covers three of the four cohorts whose raw inputs are redistributable by URL:
#   - GSE96058 / SCAN-B (GEO)
#   - METABRIC (cBioPortal datahub)
#   - GSE20685 (GEO series matrix + GPL570 annotation)
# TCGA-BRCA expression is obtained from the GDC/Xena (see data/raw/README.md);
# the derived TCGA pathway features ship in data/processed/.
# Gene sets: see download_gene_sets.sh.
#
# SHA256 values are cross-checked against the *.metadata.json files in
# data/processed/. The script fails loudly on any mismatch.

set -euo pipefail

verify() { echo "$1  $2" | shasum -a 256 -c -; }

echo "=== [1/3] GSE96058 / SCAN-B gene expression (~565MB compressed) ==="
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
