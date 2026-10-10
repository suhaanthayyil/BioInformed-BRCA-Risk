# Phase 6 External Validation QC

Generated: 2026-10-03T18:46:53-04:00
Headline model: DeepSurv

## External Validation Metrics

| cohort    | model             |    n |   events |   harrell_c |   harrell_c_ci_low |   harrell_c_ci_high | status   |
|:----------|:------------------|-----:|---------:|------------:|-------------------:|--------------------:|:---------|
| GSE20685  | Cox_PH            |  327 |       83 |       0.638 |              0.569 |               0.696 | ok       |
| GSE20685  | DeepSurv          |  327 |       83 |       0.672 |              0.607 |               0.731 | ok       |
| GSE20685  | PAM50_RORS_genefu |  327 |       83 |       0.635 |              0.566 |               0.693 | ok       |
| GSE96058  | Cox_PH            | 3273 |      336 |       0.610 |              0.579 |               0.642 | ok       |
| GSE96058  | DeepSurv          | 3273 |      336 |       0.569 |              0.535 |               0.600 | ok       |
| GSE96058  | PAM50_RORS_genefu | 3273 |      336 |       0.646 |              0.615 |               0.674 | ok       |
| METABRIC  | Cox_PH            | 1979 |     1143 |       0.564 |              0.545 |               0.583 | ok       |
| METABRIC  | DeepSurv          | 1979 |     1143 |       0.552 |              0.534 |               0.569 | ok       |
| METABRIC  | PAM50_RORS_genefu | 1979 |     1143 |       0.595 |              0.578 |               0.612 | ok       |
| TCGA-BRCA | Cox_PH            | 1064 |      150 |       0.617 |              0.564 |               0.657 | ok       |
| TCGA-BRCA | DeepSurv          | 1064 |      150 |       0.591 |              0.539 |               0.638 | ok       |
| TCGA-BRCA | PAM50_RORS_genefu | 1064 |      150 |       0.615 |              0.560 |               0.665 | ok       |

## Primary Endpoint Rows

| cohort    | subgroup   | headline_model   | comparator        |   n |   events |   delta_cindex |   ci_low |   ci_high |     p |
|:----------|:-----------|:-----------------|:------------------|----:|---------:|---------------:|---------:|----------:|------:|
| GSE96058  | tnbc       | DeepSurv         | PAM50_RORS_genefu | 143 |       26 |          0.002 |   -0.150 |     0.158 | 0.980 |
| METABRIC  | tnbc       | DeepSurv         | PAM50_RORS_genefu | 319 |      168 |          0.016 |   -0.052 |     0.090 | 0.657 |
| TCGA-BRCA | tnbc       | DeepSurv         | PAM50_RORS_genefu | 157 |       25 |          0.079 |   -0.115 |     0.253 | 0.405 |

## Random-Effects Primary Endpoint

| generated_at              | headline_model   | primary_endpoint                                       | pam50_comparator   |   threshold_delta |   threshold_p | met   |   meta_n_studies |   meta_effect |   meta_se |   meta_ci_low |   meta_ci_high |   meta_z |   meta_p |   meta_tau2 |   meta_i2 | cohorts_in_meta                       |
|:--------------------------|:-----------------|:-------------------------------------------------------|:-------------------|------------------:|--------------:|:------|-----------------:|--------------:|----------:|--------------:|---------------:|---------:|---------:|------------:|----------:|:--------------------------------------|
| 2026-10-03T18:46:53-04:00 | DeepSurv         | TNBC random-effects delta C-index vs PAM50_RORS_genefu | PAM50_RORS_genefu  |            0.0300 |        0.0500 | False |                3 |        0.0204 |    0.0310 |       -0.0403 |         0.0811 |   0.6595 |   0.5096 |      0.0000 |    0.0000 | ['GSE96058', 'METABRIC', 'TCGA-BRCA'] |
