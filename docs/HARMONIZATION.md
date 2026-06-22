# Cohort harmonization (reviewer R1.4)

How the four cohorts were harmonized into a single modeling frame, and the
residual cross-platform heterogeneity that remains. All numbers below are
reproducible from `results/reports/cohort_gene_overlap.csv`,
`results/reports/pathway_scoring_coverage.csv`,
`results/reports/harmonization_summary.json`, and
`data/processed/*.metadata.json`.

## Platforms and probe→gene mapping

| Cohort | Platform | Genes | Retained n | Probe→gene collapse |
|--------|----------|------:|-----------:|---------------------|
| TCGA-BRCA | RNA-seq, GDC log2 RSEM-normalized | (GDC) | 213 | gene-level already (HGNC symbols) |
| GSE96058 / SCAN-B | RNA-seq, GEO transformed | 30,865 | 1,483 | gene-level (symbols) |
| METABRIC | Illumina microarray, cBioPortal z-scores vs diploid | 20,603 | 1,981 | gene-level (cBioPortal collapsed) |
| GSE20685 | Affymetrix HG-U133 Plus 2.0 (GPL570) | 22,189 | 327 | probes mapped to symbols via `data/raw/GPL570.annot.gz`; multi-probe genes collapsed to the max-variance probe |

Gene symbols are upper-cased and matched to MSigDB v2024.1.Hs gene-set symbols
(`scripts/score_pathways.py`). 9,509 GSE20685 probes without a gene symbol were
dropped (recorded in `data/processed/04_gse20685.metadata.json`).

## Cross-cohort gene overlap (pathway-set coverage)

Pathway scoring requires each MSigDB set's member genes to be present. Mean gene
coverage across the 347 scored sets is high in every cohort, so the seven
harmonized features are computed from near-complete gene membership rather than
from a small lowest-common-denominator gene set:

| Cohort | Pathway sets scored | Mean gene coverage | Min set coverage |
|--------|--------------------:|-------------------:|-----------------:|
| METABRIC | 347 | 0.985 | 0.349 |
| GSE96058 | 347 | 0.964 | 0.106 |
| TCGA-BRCA | 345 | 0.947 | 0.106 |
| GSE20685 | 344 | 0.944 | 0.106 |

Per-set, per-cohort coverage is in `results/reports/pathway_scoring_coverage.csv`.

## Normalization — no cross-cohort leakage

1. **Pathway-set scores** are computed per cohort, then z-scored **within each
   cohort** (`aggregate_seven_pathways` in `src/pathways.py`:
   `z = (z - z.mean) / z.std` over that cohort's samples). The seven features are
   the mean of their constituent set z-scores. Because standardization is
   within-cohort, raw platform scale differences (e.g. METABRIC microarray
   z-scores vs RNA-seq) do not propagate across cohorts.
2. **Model preprocessing** (`SimpleImputer(median)` → `StandardScaler`) is fit on
   **TCGA-BRCA training data only** and applied unchanged (`.transform`) to every
   external cohort (`scripts/external_validation.py`). No external label or
   distribution information enters training — there is no cross-cohort leakage.

## Residual heterogeneity

Even after harmonization, platform and population differences remain. This is
quantified directly in the random-effects meta-analysis of the TNBC primary
endpoint: I² = 10.5%, τ² = 0.00045
(`data/processed/phase6_primary_endpoint.json`) — low-to-moderate heterogeneity
across the three TNBC-evaluable cohorts. Per-cohort discrimination differs
substantially (e.g. external GBSA Harrell C ranges 0.55–0.66 across cohorts;
`results/Table_1_ml_external_validation.csv`), which is consistent with genuine
cross-cohort heterogeneity and is the reason external superiority over PAM50-ROR
was not established.
