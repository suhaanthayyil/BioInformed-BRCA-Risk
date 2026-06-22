#!/bin/bash
# Download the pinned MSigDB v2024.1.Hs gene-set (.gmt) files used by
# scripts/score_pathways.py (MSIGDB_VERSION = "v2024.1.Hs").
# Run from data/raw/.  Verifies SHA256 against the values recorded in
# data/processed/msigdb_v2024_1_Hs.metadata.json.
#
# NOTE: MSigDB may require (free) registration at https://www.gsea-msigdb.org/.
# If a direct download returns HTML/login instead of a .gmt, download the files
# manually from the release directory below and place them in data/raw/gene_sets/.

set -euo pipefail

BASE="https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs"
DEST="gene_sets"
mkdir -p "$DEST"

# file<TAB>sha256
read -r -d '' MANIFEST <<'EOF' || true
h.all.v2024.1.Hs.symbols.gmt	ee2463540042078bfa3f67828e1e223bb354446d9fbb4d22845866835ba5c772
c2.cp.reactome.v2024.1.Hs.symbols.gmt	9ea1b5e656597daf423e41c5ebcaa9892bfedf3292fff768605d7b0d5e5e9703
c2.cp.kegg_medicus.v2024.1.Hs.symbols.gmt	0b4e0d7c05e8fb12f9ffcc42bc69e54a590c6d68cabb4062f17bdc37ac57674b
c2.cp.pid.v2024.1.Hs.symbols.gmt	208772254b8c13a6ee1c4a864cfa14b56657f279ec5b36e8999db87b951155e5
c5.go.bp.v2024.1.Hs.symbols.gmt	71df041c352cc810c11847db328bb94925c5448a009f310d224f2e8baa60b792
EOF

while IFS=$'\t' read -r fname sha; do
  [ -z "$fname" ] && continue
  out="$DEST/$fname"
  if [ -f "$out" ]; then
    echo "exists: $out"
  else
    echo "=== downloading $fname ==="
    curl -fL -o "$out" "$BASE/$fname"
  fi
  echo "$sha  $out" | shasum -a 256 -c -
done <<< "$MANIFEST"

echo "=== All gene-set files present and verified ==="
