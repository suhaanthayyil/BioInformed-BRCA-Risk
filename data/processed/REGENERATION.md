# Regenerating processed data from raw inputs

Two reproducibility tiers are supported.

## Tier 1 — reproduce the headline results without raw data (recommended for review)

The derived intermediates the modeling pipeline actually consumes are committed
to the repository, so `scripts/train_ml_zoo.py`, `scripts/external_validation.py`,
and all Phase-3 analyses run on a fresh clone with **no multi-GB raw download**:

| Committed intermediate | Consumed by |
|------------------------|-------------|
| `data/processed/04_pathway_features.parquet` | train, external validation, stability, ablation |
| `data/processed/pathway_scores_all.parquet` | pathway-scoring sensitivity |
| `data/processed/baselines_pam50.parquet` | external validation (PAM50-ROR comparator) |
| `data/processed/unified_cohorts.duckdb` | harmonized clinical/survival for all cohorts |
| `data/processed/ml_model_zoo.metadata.json` | locked feature list + selected configs |
| `data/processed/02_tcga_feature_matrix.csv` | canonical TCGA sample-ID ordering |
| `models/*.pkl`, `models/*.pt` | the locked trained weights |

Verify with `python scripts/check_submission_consistency.py` and `make smoke`.

## Tier 2 — rebuild everything from raw expression

The large raw expression matrices (`01_tcga_expression_normalized.csv`,
`02_gse96058_expression.parquet`, `03_metabric_expression.parquet`,
`04_gse20685_expression.parquet`) are **not redistributed** here (size + source
terms; see `docs/DATA_USE.md`). To rebuild from scratch:

1. **Download raw inputs** (with checksum verification):
   ```bash
   cd data/raw
   bash download_data.sh         # GSE96058, METABRIC, GSE20685
   bash download_gene_sets.sh    # MSigDB v2024.1.Hs .gmt files
   # TCGA-BRCA expression: obtain the log2-normalized RSEM matrix from
   # GDC/Xena (see data/raw/README.md) -> data/01_tcga_expression_normalized.csv
   ```

2. **Parse raw → expression parquets.** This step currently lives in the
   exploratory notebooks (run top-to-bottom):
   - `notebooks/01_data_preprocessing.ipynb` → `02_gse96058_expression.parquet`,
     `03_metabric_expression.parquet`, `04_gse20685_expression.parquet`, and the
     TCGA sample-ID alignment behind `02_tcga_feature_matrix.csv`.
   These notebooks are provided for transparency; porting them into a single
   headless ingestion script is tracked as follow-up work. Until then, Tier 1
   is the supported review path.

3. **Rebuild the pipeline** from the parsed expression:
   ```bash
   make all     # harmonize -> score_pathways -> baselines -> train -> external_validation -> phase3
   ```

`*.metadata.json` files in this directory record the exact source URL, SHA256,
and shape of every raw input, so a rebuild is fully traceable.
