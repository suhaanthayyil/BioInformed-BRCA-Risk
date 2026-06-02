# Phase 3 Clinical Baseline QC

Generated: 2026-05-16T13:11:08-04:00

## Score Inventory

| cohort    | baseline                      |   n_scored |
|:----------|:------------------------------|-----------:|
| GSE20685  | PAM50_ROR_official            |        327 |
| GSE96058  | Oncotype_DX_21_gene_surrogate |       1251 |
| GSE96058  | PAM50_ROR_official            |       1483 |
| GSE96058  | PAM50_ROR_surrogate           |       1480 |
| METABRIC  | Oncotype_DX_21_gene_surrogate |       1825 |
| METABRIC  | PAM50_ROR_official            |       1980 |
| METABRIC  | PAM50_ROR_surrogate           |       1948 |
| TCGA-BRCA | Oncotype_DX_21_gene_surrogate |        118 |
| TCGA-BRCA | PAM50_ROR_official            |        213 |
| TCGA-BRCA | PAM50_ROR_surrogate           |        211 |

## Unavailable / Surrogate Notes

| cohort    | baseline                      | status                          | note                                                                                                                                                                                                                                                                                                                                               |
|:----------|:------------------------------|:--------------------------------|:---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| GSE20685  | PAM50_ROR_surrogate           | unavailable                     | PAM50-ROR surrogate, not official ROR-S/ROR-P/ROR-PT; published intrinsic subtype risk mapping, proliferation modifier (5 genes). Official genefu PAM50-ROR is the primary comparator; this surrogate is retained only for historical comparison.                                                                                                                                    |
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
| GSE20685  | PAM50_ROR_official            |  327 |       83 |       0.635 |              0.566 |               0.693 |   0.600 | ok             |
| GSE20685  | PAM50_ROR_surrogate           |    0 |        0 |     nan     |            nan     |             nan     | nan     | too_few_events |
| GSE96058  | Oncotype_DX_21_gene_surrogate | 1251 |      242 |       0.567 |              0.530 |               0.602 |   0.577 | ok             |
| GSE96058  | PAM50_ROR_official            | 1483 |      322 |       0.633 |              0.603 |               0.662 |   0.633 | ok             |
| GSE96058  | PAM50_ROR_surrogate           | 1480 |      321 |       0.638 |              0.604 |               0.670 |   0.638 | ok             |
| METABRIC  | Oncotype_DX_21_gene_surrogate | 1507 |      882 |       0.587 |              0.566 |               0.609 |   0.579 | ok             |
| METABRIC  | PAM50_ROR_official            | 1980 |     1143 |       0.595 |              0.577 |               0.613 |   0.583 | ok             |
| METABRIC  | PAM50_ROR_surrogate           | 1948 |     1122 |       0.588 |              0.570 |               0.607 |   0.574 | ok             |
| TCGA-BRCA | Oncotype_DX_21_gene_surrogate |  118 |       62 |       0.578 |              0.483 |               0.669 |   0.571 | ok             |
| TCGA-BRCA | PAM50_ROR_official            |  213 |      132 |       0.504 |              0.445 |               0.556 |   0.507 | ok             |
| TCGA-BRCA | PAM50_ROR_surrogate           |  211 |      130 |       0.552 |              0.490 |               0.611 |   0.553 | ok             |

## TNBC OS Metrics

| cohort    | baseline                      |   n |   events |   harrell_c |   harrell_c_ci_low |   harrell_c_ci_high | status         |
|:----------|:------------------------------|----:|---------:|------------:|-------------------:|--------------------:|:---------------|
| GSE20685  | Oncotype_DX_21_gene_surrogate |   0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| GSE20685  | PAM50_ROR_official            |   0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| GSE20685  | PAM50_ROR_surrogate           |   0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| GSE96058  | Oncotype_DX_21_gene_surrogate |   0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| GSE96058  | PAM50_ROR_official            |  55 |       27 |       0.576 |              0.468 |               0.686 | ok             |
| GSE96058  | PAM50_ROR_surrogate           |  55 |       27 |       0.584 |              0.467 |               0.691 | ok             |
| METABRIC  | Oncotype_DX_21_gene_surrogate |   0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| METABRIC  | PAM50_ROR_official            | 320 |      168 |       0.525 |              0.481 |               0.571 | ok             |
| METABRIC  | PAM50_ROR_surrogate           | 312 |      163 |       0.581 |              0.535 |               0.625 | ok             |
| TCGA-BRCA | Oncotype_DX_21_gene_surrogate |   0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| TCGA-BRCA | PAM50_ROR_official            |  33 |       19 |       0.502 |              0.371 |               0.643 | ok             |
| TCGA-BRCA | PAM50_ROR_surrogate           |  33 |       19 |       0.490 |              0.353 |               0.624 | ok             |

PAM50 values are marked as `PAM50_ROR_surrogate` unless official genefu scoring is available. MammaPrint, EndoPredict, and BCI are not fabricated when exact public formulae are absent.
