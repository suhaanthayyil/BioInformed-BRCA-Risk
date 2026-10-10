# Phase 3 Clinical Baseline QC

Generated: 2026-10-03T17:36:42-04:00

## Score Inventory

| cohort    | baseline                      |   n_scored |
|:----------|:------------------------------|-----------:|
| GSE20685  | PAM50_RORS_genefu             |        327 |
| GSE96058  | Oncotype_DX_21_gene_surrogate |       2832 |
| GSE96058  | PAM50_RORS_genefu             |       3273 |
| GSE96058  | PAM50_ROR_surrogate           |       3239 |
| METABRIC  | Oncotype_DX_21_gene_surrogate |       1506 |
| METABRIC  | PAM50_RORS_genefu             |       1980 |
| METABRIC  | PAM50_ROR_surrogate           |       1948 |
| TCGA-BRCA | Oncotype_DX_21_gene_surrogate |        780 |
| TCGA-BRCA | PAM50_RORS_genefu             |       1064 |
| TCGA-BRCA | PAM50_ROR_surrogate           |        820 |

## Unavailable / Surrogate Notes

| cohort    | baseline                      | status                          | note                                                                                                                                                                                                                                                                                                                                               |
|:----------|:------------------------------|:--------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| GSE20685  | PAM50_ROR_surrogate           | unavailable                     | PAM50-ROR surrogate, not genefu ROR-S/ROR-P/ROR-PT; published intrinsic subtype risk mapping, proliferation modifier (5 genes). genefu PAM50 ROR-S (research implementation) is the primary comparator; this surrogate is retained only for historical comparison.                                                                                 |
| GSE20685  | Oncotype_DX_21_gene_surrogate | not_evaluable_missing_er_status | Paik 2004 group-weight surrogate on within-cohort z-scored expression; genes present ER=4, proliferation=5, HER2=2, invasion=2. Commercial reference normalization is unavailable. Restricted to ER-positive samples for performance evaluation. No ER status is available for this cohort, so the ER-positive Oncotype analysis is not evaluable. |
| GSE20685  | MammaPrint_gene70             | unavailable                     | Exact MammaPrint scoring was not computed in Python because the public repo does not contain the validated 70-gene coefficients/centroid. Use genefu::gene70 if the R package installation succeeds; otherwise report as not reproducible from public inputs rather than fabricating a score.                                                      |
| GSE20685  | EndoPredict                   | unavailable                     | EndoPredict exact public formula/coefficients were not available in the repository; score not fabricated.                                                                                                                                                                                                                                          |
| GSE20685  | Breast_Cancer_Index           | unavailable                     | Breast_Cancer_Index exact public formula/coefficients were not available in the repository; score not fabricated.                                                                                                                                                                                                                                  |
| GSE96058  | MammaPrint_gene70             | unavailable                     | Exact MammaPrint scoring was not computed in Python because the public repo does not contain the validated 70-gene coefficients/centroid. Use genefu::gene70 if the R package installation succeeds; otherwise report as not reproducible from public inputs rather than fabricating a score.                                                      |
| GSE96058  | EndoPredict                   | unavailable                     | EndoPredict exact public formula/coefficients were not available in the repository; score not fabricated.                                                                                                                                                                                                                                          |
| GSE96058  | Breast_Cancer_Index           | unavailable                     | Breast_Cancer_Index exact public formula/coefficients were not available in the repository; score not fabricated.                                                                                                                                                                                                                                  |
| METABRIC  | MammaPrint_gene70             | unavailable                     | Exact MammaPrint scoring was not computed in Python because the public repo does not contain the validated 70-gene coefficients/centroid. Use genefu::gene70 if the R package installation succeeds; otherwise report as not reproducible from public inputs rather than fabricating a score.                                                      |
| METABRIC  | EndoPredict                   | unavailable                     | EndoPredict exact public formula/coefficients were not available in the repository; score not fabricated.                                                                                                                                                                                                                                          |
| METABRIC  | Breast_Cancer_Index           | unavailable                     | Breast_Cancer_Index exact public formula/coefficients were not available in the repository; score not fabricated.                                                                                                                                                                                                                                  |
| TCGA-BRCA | MammaPrint_gene70             | unavailable                     | Exact MammaPrint scoring was not computed in Python because the public repo does not contain the validated 70-gene coefficients/centroid. Use genefu::gene70 if the R package installation succeeds; otherwise report as not reproducible from public inputs rather than fabricating a score.                                                      |
| TCGA-BRCA | EndoPredict                   | unavailable                     | EndoPredict exact public formula/coefficients were not available in the repository; score not fabricated.                                                                                                                                                                                                                                          |
| TCGA-BRCA | Breast_Cancer_Index           | unavailable                     | Breast_Cancer_Index exact public formula/coefficients were not available in the repository; score not fabricated.                                                                                                                                                                                                                                  |

## Overall OS Metrics

| cohort    | baseline                      |    n |   events |   harrell_c |   harrell_c_ci_low |   harrell_c_ci_high |   uno_c | status         |
|:----------|:------------------------------|-----:|---------:|------------:|-------------------:|--------------------:|--------:|:---------------|
| GSE20685  | Oncotype_DX_21_gene_surrogate |    0 |        0 |     nan     |            nan     |             nan     | nan     | too_few_events |
| GSE20685  | PAM50_RORS_genefu             |  327 |       83 |       0.635 |              0.566 |               0.693 |   0.600 | ok             |
| GSE20685  | PAM50_ROR_surrogate           |    0 |        0 |     nan     |            nan     |             nan     | nan     | too_few_events |
| GSE96058  | Oncotype_DX_21_gene_surrogate | 2832 |      255 |       0.567 |              0.529 |               0.604 |   0.576 | ok             |
| GSE96058  | PAM50_RORS_genefu             | 3273 |      336 |       0.646 |              0.615 |               0.674 |   0.628 | ok             |
| GSE96058  | PAM50_ROR_surrogate           | 3239 |      335 |       0.650 |              0.616 |               0.680 |   0.640 | ok             |
| METABRIC  | Oncotype_DX_21_gene_surrogate | 1506 |      881 |       0.587 |              0.565 |               0.609 |   0.579 | ok             |
| METABRIC  | PAM50_RORS_genefu             | 1980 |     1143 |       0.595 |              0.577 |               0.613 |   0.583 | ok             |
| METABRIC  | PAM50_ROR_surrogate           | 1948 |     1122 |       0.589 |              0.570 |               0.608 |   0.574 | ok             |
| TCGA-BRCA | Oncotype_DX_21_gene_surrogate |  780 |       99 |       0.574 |              0.483 |               0.644 |   0.533 | ok             |
| TCGA-BRCA | PAM50_RORS_genefu             | 1064 |      150 |       0.615 |              0.560 |               0.665 |   0.577 | ok             |
| TCGA-BRCA | PAM50_ROR_surrogate           |  820 |      125 |       0.569 |              0.507 |               0.631 |   0.528 | ok             |

## TNBC OS Metrics

| cohort    | baseline                      |   n |   events |   harrell_c |   harrell_c_ci_low |   harrell_c_ci_high | status         |
|:----------|:------------------------------|----:|---------:|------------:|-------------------:|--------------------:|:---------------|
| GSE20685  | Oncotype_DX_21_gene_surrogate |   0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| GSE20685  | PAM50_RORS_genefu             |   0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| GSE20685  | PAM50_ROR_surrogate           |   0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| GSE96058  | Oncotype_DX_21_gene_surrogate |   0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| GSE96058  | PAM50_RORS_genefu             | 143 |       26 |       0.577 |              0.487 |               0.667 | ok             |
| GSE96058  | PAM50_ROR_surrogate           | 139 |       26 |       0.585 |              0.463 |               0.700 | ok             |
| METABRIC  | Oncotype_DX_21_gene_surrogate |   0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| METABRIC  | PAM50_RORS_genefu             | 320 |      168 |       0.525 |              0.481 |               0.571 | ok             |
| METABRIC  | PAM50_ROR_surrogate           | 312 |      163 |       0.579 |              0.533 |               0.622 | ok             |
| TCGA-BRCA | Oncotype_DX_21_gene_surrogate |   0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| TCGA-BRCA | PAM50_RORS_genefu             | 157 |       25 |       0.440 |              0.304 |               0.566 | ok             |
| TCGA-BRCA | PAM50_ROR_surrogate           | 113 |       16 |       0.421 |              0.283 |               0.585 | ok             |

PAM50 values are marked as `PAM50_ROR_surrogate` unless genefu scoring is available. MammaPrint, EndoPredict, and BCI are not fabricated when exact public formulae are absent.
