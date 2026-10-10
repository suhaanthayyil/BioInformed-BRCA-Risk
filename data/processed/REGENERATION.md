# Regenerating processed data from raw inputs

Two reproducibility tiers are supported.

## Tier 1 - reproduce the results without raw data (recommended for review)

The intermediates the modeling pipeline consumes are committed, so
`scripts/train_ml_zoo.py`, `scripts/external_validation.py` and every
downstream analysis run on a fresh clone with no raw download:

| Committed intermediate | Consumed by |
|------------------------|-------------|
| `data/processed/04_pathway_features.parquet` | train, external validation, stability, ablation |
| `data/processed/pathway_scores_all.parquet` | pathway-scoring sensitivity |
| `data/processed/baselines_pam50.parquet` | genefu PAM50 ROR-S comparator |
| `data/processed/baselines_oncotype.parquet` | Oncotype DX surrogate comparison (Table S18) |
| `data/processed/unified_cohorts.duckdb` | harmonized clinical/survival for all cohorts |
| `data/processed/ingest_flow.json` | cohort inclusion counts (Figure 1) |
| `data/processed/ml_model_zoo.metadata.json` | feature list + selected configs |
| `models/*.pkl`, `models/*.pt` | trained weights; `models/locked_v1/` holds the 2026-05-16 weights |

Run `make all`, then `make test` and `make consistency`.

## Tier 2 - rebuild everything from raw expression

Raw expression matrices are not redistributed (size and source terms; see
`docs/DATA_USE.md`). To rebuild from scratch (needs R with genefu):

```bash
cd data/raw
bash download_data.sh         # TCGA-BRCA (Xena + TCGA-CDR), GSE96058, METABRIC, GSE20685
bash download_gene_sets.sh    # MSigDB v2024.1.Hs .gmt files
cd ../..
make all-from-raw             # ingest_raw -> harmonize -> score_pathways -> baselines -> make all
```

`scripts/ingest_raw.py` applies the inclusion criteria in
`docs/COHORT_NOTES.md` and records SHA256 of every raw input in
`data/processed/ingest_flow.json`. The GSE20685 expression and clinical
parquets (`04_gse20685_*.parquet`) were parsed from the GEO series matrix and
GPL570 annotation listed in `04_gse20685.metadata.json` and are shipped.
