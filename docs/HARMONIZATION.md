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
| TCGA-BRCA | RNA-seq, UCSC Xena HiSeqV2 log2(RSEM+1) | 20,530 | 1,064 | gene-level already (HGNC symbols) |
| GSE96058 / SCAN-B | RNA-seq, GEO transformed | 30,865 | 3,273 | gene-level (symbols) |
| METABRIC | Illumina microarray, cBioPortal z-scores vs diploid | 20,603 | 1,980 analysed (of 2,509 clinical records; 1,981 records have OS, 1,980 of them have expression) | gene-level (cBioPortal collapsed) |
| GSE20685 | Affymetrix HG-U133 Plus 2.0 (GPL570) | 22,189 | 327 | probes mapped to symbols via `data/raw/GPL570.annot.gz`; multi-probe genes collapsed by averaging probes |

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
| METABRIC | 347 | 0.984 | 0.349 |
| GSE96058 | 347 | 0.964 | 0.106 |
| TCGA-BRCA | 345 | 0.947 | 0.106 |
| GSE20685 | 344 | 0.944 | 0.106 |

Per-set, per-cohort coverage is in `results/reports/pathway_scoring_coverage.csv`.

## Normalization

1. **Gene-level standardization within each cohort.** Before scoring, every
   cohort's expression is z-scored per gene over that cohort's own samples
   (`zscore_by_gene` in `src/pathways.py`, called by `scripts/score_pathways.py`).
   This includes the external cohorts. It is unsupervised (no outcome or label is
   used) and is the only step that uses an external cohort's own data; it puts
   RNA-seq and microarray platforms on a common relative scale.
2. **Pathway-set scores** are computed per sample from within-sample ranks, then
   z-scored within each cohort and averaged into the seven features
   (`aggregate_seven_pathways`).
3. **Model preprocessing** (`SimpleImputer(median)` → `StandardScaler`) and every
   model parameter are fit on **TCGA-BRCA only** and applied unchanged
   (`.transform` / `.predict`) to the external cohorts. No external outcome
   information enters training, tuning or model selection; "applied without
   refitting" refers to these fitted components.

## Residual heterogeneity

Platform and population differences remain after harmonization. For the TNBC
primary endpoint the random-effects heterogeneity is I² = 0%
(`data/processed/phase6_primary_endpoint.json`). Overall discrimination of the
selected model differs across cohorts (Harrell C 0.552 to 0.672;
`results/Table_2_external_validation.csv`), consistent with genuine platform and
population heterogeneity.
