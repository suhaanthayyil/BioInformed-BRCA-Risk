# Phase 6 External Validation QC

Generated: 2026-05-16T13:17:28-04:00
Headline model: Gradient_Boosted_Survival

## External Validation Metrics

| cohort    | model                     |    n |   events |   harrell_c |   harrell_c_ci_low |   harrell_c_ci_high | status   |
|:----------|:--------------------------|-----:|---------:|------------:|-------------------:|--------------------:|:---------|
| GSE20685  | Cox_PH                    |  327 |       83 |       0.534 |              0.472 |               0.597 | ok       |
| GSE20685  | Gradient_Boosted_Survival |  327 |       83 |       0.565 |              0.515 |               0.617 | ok       |
| GSE20685  | PAM50_ROR_official        |  327 |       83 |       0.635 |              0.566 |               0.693 | ok       |
| GSE96058  | Cox_PH                    | 1483 |      322 |       0.643 |              0.613 |               0.671 | ok       |
| GSE96058  | Gradient_Boosted_Survival | 1483 |      322 |       0.663 |              0.633 |               0.691 | ok       |
| GSE96058  | PAM50_ROR_official        | 1483 |      322 |       0.633 |              0.603 |               0.662 | ok       |
| METABRIC  | Cox_PH                    | 1979 |     1143 |       0.563 |              0.545 |               0.581 | ok       |
| METABRIC  | Gradient_Boosted_Survival | 1979 |     1143 |       0.555 |              0.540 |               0.570 | ok       |
| METABRIC  | PAM50_ROR_official        | 1979 |     1143 |       0.595 |              0.578 |               0.612 | ok       |
| TCGA-BRCA | Cox_PH                    |  213 |      132 |       0.634 |              0.582 |               0.684 | ok       |
| TCGA-BRCA | Gradient_Boosted_Survival |  213 |      132 |       0.673 |              0.628 |               0.719 | ok       |
| TCGA-BRCA | PAM50_ROR_official        |  213 |      132 |       0.504 |              0.444 |               0.557 | ok       |

## Primary Endpoint Rows

| cohort    | subgroup   | headline_model            | comparator         |   n |   events |   delta_cindex |   ci_low |   ci_high |     p |
|:----------|:-----------|:--------------------------|:-------------------|----:|---------:|---------------:|---------:|----------:|------:|
| GSE96058  | tnbc       | Gradient_Boosted_Survival | PAM50_ROR_official |  55 |       27 |          0.030 |   -0.131 |     0.195 | 0.729 |
| METABRIC  | tnbc       | Gradient_Boosted_Survival | PAM50_ROR_official | 320 |      168 |         -0.010 |   -0.070 |     0.051 | 0.744 |
| TCGA-BRCA | tnbc       | Gradient_Boosted_Survival | PAM50_ROR_official |  33 |       19 |          0.111 |   -0.027 |     0.272 | 0.151 |

## Random-Effects Primary Endpoint

| generated_at              | headline_model            | primary_endpoint                                        | pam50_comparator   |   threshold_delta |   threshold_p | met   |   meta_n_studies |   meta_effect |   meta_se |   meta_ci_low |   meta_ci_high |   meta_z |   meta_p |   meta_tau2 |   meta_i2 | cohorts_in_meta                       |
|:--------------------------|:--------------------------|:--------------------------------------------------------|:-------------------|------------------:|--------------:|:------|-----------------:|--------------:|----------:|--------------:|---------------:|---------:|---------:|------------:|----------:|:--------------------------------------|
| 2026-05-16T13:17:28-04:00 | Gradient_Boosted_Survival | TNBC random-effects delta C-index vs PAM50_ROR_official | PAM50_ROR_official |            0.0300 |        0.0500 | False |                3 |        0.0144 |    0.0314 |       -0.0472 |         0.0760 |   0.4585 |   0.6466 |      0.0004 |    0.1053 | ['GSE96058', 'METABRIC', 'TCGA-BRCA'] |
