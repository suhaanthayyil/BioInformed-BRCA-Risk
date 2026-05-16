# Phase 6 External Validation QC

Generated: 2026-05-16T12:45:44-04:00
Headline model: Gradient_Boosted_Survival

## External Validation Metrics

| cohort    | model                     |    n |   events |   harrell_c |   harrell_c_ci_low |   harrell_c_ci_high | status         |
|:----------|:--------------------------|-----:|---------:|------------:|-------------------:|--------------------:|:---------------|
| GSE20685  | Cox_PH                    |  327 |       83 |       0.534 |              0.472 |               0.597 | ok             |
| GSE20685  | Gradient_Boosted_Survival |  327 |       83 |       0.565 |              0.515 |               0.617 | ok             |
| GSE20685  | PAM50_ROR_surrogate       |    0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| GSE96058  | Cox_PH                    | 1483 |      322 |       0.643 |              0.613 |               0.671 | ok             |
| GSE96058  | Gradient_Boosted_Survival | 1483 |      322 |       0.663 |              0.633 |               0.691 | ok             |
| GSE96058  | PAM50_ROR_surrogate       | 1480 |      321 |       0.638 |              0.604 |               0.670 | ok             |
| METABRIC  | Cox_PH                    | 1979 |     1143 |       0.563 |              0.545 |               0.581 | ok             |
| METABRIC  | Gradient_Boosted_Survival | 1979 |     1143 |       0.555 |              0.540 |               0.570 | ok             |
| METABRIC  | PAM50_ROR_surrogate       | 1948 |     1122 |       0.588 |              0.568 |               0.607 | ok             |
| TCGA-BRCA | Cox_PH                    |  213 |      132 |       0.634 |              0.582 |               0.684 | ok             |
| TCGA-BRCA | Gradient_Boosted_Survival |  213 |      132 |       0.673 |              0.628 |               0.719 | ok             |
| TCGA-BRCA | PAM50_ROR_surrogate       |  211 |      130 |       0.552 |              0.490 |               0.608 | ok             |

## Primary Endpoint Rows

| cohort    | subgroup   | headline_model            | comparator          |   n |   events |   delta_cindex |   ci_low |   ci_high |     p |
|:----------|:-----------|:--------------------------|:--------------------|----:|---------:|---------------:|---------:|----------:|------:|
| GSE96058  | tnbc       | Gradient_Boosted_Survival | PAM50_ROR_surrogate |  55 |       27 |          0.022 |   -0.135 |     0.188 | 0.801 |
| METABRIC  | tnbc       | Gradient_Boosted_Survival | PAM50_ROR_surrogate | 312 |      163 |         -0.061 |   -0.125 |     0.001 | 0.055 |
| TCGA-BRCA | tnbc       | Gradient_Boosted_Survival | PAM50_ROR_surrogate |  33 |       19 |          0.123 |   -0.051 |     0.322 | 0.201 |

## Random-Effects Primary Endpoint

| generated_at              | headline_model            | primary_endpoint                                         |   threshold_delta |   threshold_p | met   |   meta_n_studies |   meta_effect |   meta_se |   meta_ci_low |   meta_ci_high |   meta_z |   meta_p |   meta_tau2 |   meta_i2 | cohorts_in_meta                       |
|:--------------------------|:--------------------------|:---------------------------------------------------------|------------------:|--------------:|:------|-----------------:|--------------:|----------:|--------------:|---------------:|---------:|---------:|------------:|----------:|:--------------------------------------|
| 2026-05-16T12:45:44-04:00 | Gradient_Boosted_Survival | TNBC random-effects delta C-index vs PAM50_ROR_surrogate |            0.0300 |        0.0500 | False |                3 |       -0.0005 |    0.0537 |       -0.1057 |         0.1047 |  -0.0090 |   0.9928 |      0.0044 |    0.4895 | ['GSE96058', 'METABRIC', 'TCGA-BRCA'] |
