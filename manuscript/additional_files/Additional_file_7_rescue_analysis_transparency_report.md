# Additional File 7. Rescue Analysis Transparency Report

This analysis informed limitations only and is not part of the locked primary endpoint.

# Rescue Analysis QC

Generated: 2026-05-16T15:48:56-04:00

This is an exploratory transparent rescue run. The locked Phase 6 primary endpoint remains unchanged.

## METABRIC Expression Audit

| current_processed_expression                    | raw_expression_file_present   | zscore_expression_file_present   | raw_shape_header_check   |   raw_first_rows_median |   raw_first_rows_sd | zscore_shape_header_check   |   zscore_first_rows_median |   zscore_first_rows_sd |
|:------------------------------------------------|:------------------------------|:---------------------------------|:-------------------------|------------------------:|--------------------:|:----------------------------|---------------------------:|-----------------------:|
| cBioPortal z-scores relative to diploid samples | True                          | True                             | [5, 1982]                |                 6.88667 |             1.38935 | [5, 1982]                   |                    -0.1198 |                      1 |

## Top TCGA-Selected Candidates

| feature_set                 | model                     |   n_features |   overall_cv_cindex |   tnbc_cv_cindex |   composite_cv_cindex |   valid_tnbc_folds |
|:----------------------------|:--------------------------|-------------:|--------------------:|-----------------:|----------------------:|-------------------:|
| component7_interactions     | DeepSurv                  |           43 |               0.635 |            0.648 |                 0.641 |                  4 |
| ensemble_top3               | Rank_Average_Ensemble     |            0 |               0.645 |            0.628 |                 0.636 |                  1 |
| seven_interactions          | Elastic_Net_Cox           |           14 |               0.623 |            0.617 |                 0.620 |                  4 |
| seven_interactions          | Random_Survival_Forest    |           14 |               0.603 |            0.634 |                 0.619 |                  4 |
| seven_interactions          | Cox_PH                    |           14 |               0.623 |            0.611 |                 0.617 |                  4 |
| seven_interactions          | Elastic_Net_Cox           |           14 |               0.619 |            0.611 |                 0.615 |                  4 |
| seven_interactions          | Elastic_Net_Cox           |           14 |               0.619 |            0.611 |                 0.615 |                  4 |
| seven_interactions          | Random_Survival_Forest    |           14 |               0.610 |            0.612 |                 0.611 |                  4 |
| seven_interactions          | DeepSurv                  |           14 |               0.628 |            0.589 |                 0.608 |                  4 |
| seven_interactions          | DeepSurv                  |           14 |               0.631 |            0.578 |                 0.604 |                  4 |
| seven_interactions          | Random_Survival_Forest    |           14 |               0.611 |            0.597 |                 0.604 |                  4 |
| seven_interactions          | Elastic_Net_Cox           |           14 |               0.623 |            0.583 |                 0.603 |                  4 |
| seven_interactions          | Gradient_Boosted_Survival |           14 |               0.592 |            0.607 |                 0.600 |                  4 |
| component7_interactions     | DeepSurv                  |           43 |               0.642 |            0.548 |                 0.595 |                  4 |
| seven_receptor_interactions | Gradient_Boosted_Survival |           23 |               0.602 |            0.581 |                 0.592 |                  4 |

## Rescue Headline

Headline rescue model: `Rescue_DeepSurv` selected by TCGA-only composite CV.

## External TNBC Rows vs Official PAM50

| cohort    | subgroup   | headline_model   | comparator         |   n |   events |   delta_cindex |   ci_low |   ci_high |       p |
|:----------|:-----------|:-----------------|:-------------------|----:|---------:|---------------:|---------:|----------:|--------:|
| GSE20685  | tnbc       | Rescue_DeepSurv  | PAM50_ROR_official |   0 |        0 |        nan     |  nan     |   nan     | nan     |
| GSE96058  | tnbc       | Rescue_DeepSurv  | PAM50_ROR_official |  55 |       27 |         -0.039 |   -0.223 |     0.121 |   0.662 |
| METABRIC  | tnbc       | Rescue_DeepSurv  | PAM50_ROR_official | 320 |      168 |          0.009 |   -0.062 |     0.078 |   0.805 |
| TCGA-BRCA | tnbc       | Rescue_DeepSurv  | PAM50_ROR_official |  33 |       19 |          0.424 |    0.273 |     0.591 |   0.000 |

## Random-Effects Meta-Analysis

| label                      | cohorts                               | met   |   meta_n_studies |   meta_effect |   meta_se |   meta_ci_low |   meta_ci_high |   meta_z |   meta_p |   meta_tau2 |   meta_i2 |
|:---------------------------|:--------------------------------------|:------|-----------------:|--------------:|----------:|--------------:|---------------:|---------:|---------:|------------:|----------:|
| all_protocol_tnbc_cohorts  | ['GSE96058', 'METABRIC', 'TCGA-BRCA'] | False |                3 |        0.1284 |    0.1324 |       -0.1312 |         0.3880 |   0.9697 |   0.3322 |      0.0476 |    0.9152 |
| external_only_tnbc_cohorts | ['GSE96058', 'METABRIC']              | False |                2 |        0.0023 |    0.0330 |       -0.0624 |         0.0669 |   0.0688 |   0.9452 |      0.0000 |    0.0000 |

## External Validation Metrics for Headline and PAM50

| cohort    | model              | subgroup   |    n |   events |   harrell_c |   harrell_c_ci_low |   harrell_c_ci_high | status         |
|:----------|:-------------------|:-----------|-----:|---------:|------------:|-------------------:|--------------------:|:---------------|
| GSE20685  | Rescue_DeepSurv    | overall    |  327 |       83 |       0.585 |              0.521 |               0.645 | ok             |
| GSE20685  | Rescue_DeepSurv    | tnbc       |    0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| GSE20685  | PAM50_ROR_official | overall    |  327 |       83 |       0.635 |              0.566 |               0.693 | ok             |
| GSE20685  | PAM50_ROR_official | tnbc       |    0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| GSE96058  | Rescue_DeepSurv    | overall    | 1483 |      322 |       0.598 |              0.567 |               0.631 | ok             |
| GSE96058  | Rescue_DeepSurv    | tnbc       |   55 |       27 |       0.537 |              0.423 |               0.642 | ok             |
| GSE96058  | PAM50_ROR_official | overall    | 1483 |      322 |       0.633 |              0.603 |               0.662 | ok             |
| GSE96058  | PAM50_ROR_official | tnbc       |   55 |       27 |       0.576 |              0.468 |               0.686 | ok             |
| METABRIC  | Rescue_DeepSurv    | overall    | 1979 |     1143 |       0.582 |              0.563 |               0.599 | ok             |
| METABRIC  | Rescue_DeepSurv    | tnbc       |  319 |      168 |       0.534 |              0.485 |               0.580 | ok             |
| METABRIC  | PAM50_ROR_official | overall    | 1979 |     1143 |       0.595 |              0.578 |               0.612 | ok             |
| METABRIC  | PAM50_ROR_official | tnbc       |  319 |      168 |       0.525 |              0.481 |               0.569 | ok             |
| TCGA-BRCA | Rescue_DeepSurv    | overall    |  213 |      132 |       0.926 |              0.904 |               0.946 | ok             |
| TCGA-BRCA | Rescue_DeepSurv    | tnbc       |   33 |       19 |       0.927 |              0.868 |               0.978 | ok             |
| TCGA-BRCA | PAM50_ROR_official | overall    |  213 |      132 |       0.504 |              0.444 |               0.557 | ok             |
| TCGA-BRCA | PAM50_ROR_official | tnbc       |   33 |       19 |       0.502 |              0.355 |               0.655 | ok             |
