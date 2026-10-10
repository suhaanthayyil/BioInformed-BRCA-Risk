# Rescue Analysis QC

Generated: 2026-10-01T23:00:56-04:00

This is an exploratory post-hoc screen. It does not replace the registered primary endpoint.

## METABRIC Expression Audit

| current_processed_expression                    | raw_expression_file_present   | zscore_expression_file_present   | zscore_shape_header_check   |   zscore_first_rows_median |   zscore_first_rows_sd |
|:------------------------------------------------|:------------------------------|:---------------------------------|:----------------------------|---------------------------:|-----------------------:|
| cBioPortal z-scores relative to diploid samples | False                         | True                             | [5, 1982]                   |                    -0.1198 |                      1 |

## Top TCGA-Selected Candidates

| feature_set        | model                 |   n_features |   overall_cv_cindex |   tnbc_cv_cindex |   composite_cv_cindex |   valid_tnbc_folds |
|:-------------------|:----------------------|-------------:|--------------------:|-----------------:|----------------------:|-------------------:|
| seven              | Cox_PH                |            9 |               0.735 |            0.831 |                 0.783 |                  5 |
| seven              | Elastic_Net_Cox       |            9 |               0.734 |            0.831 |                 0.782 |                  5 |
| ensemble_top3      | Rank_Average_Ensemble |            0 |               0.744 |            0.820 |                 0.782 |                  1 |
| hallmark50         | Elastic_Net_Cox       |           52 |               0.743 |            0.820 |                 0.782 |                  5 |
| seven              | Cox_PH                |            9 |               0.734 |            0.828 |                 0.781 |                  5 |
| seven              | Cox_PH                |            9 |               0.744 |            0.811 |                 0.778 |                  5 |
| seven              | Elastic_Net_Cox       |            9 |               0.734 |            0.819 |                 0.776 |                  5 |
| seven              | Cox_PH                |            9 |               0.734 |            0.818 |                 0.776 |                  5 |
| hallmark50         | Elastic_Net_Cox       |           52 |               0.730 |            0.821 |                 0.775 |                  5 |
| seven              | Elastic_Net_Cox       |            9 |               0.732 |            0.818 |                 0.775 |                  5 |
| seven              | Cox_PH                |            9 |               0.750 |            0.798 |                 0.774 |                  5 |
| seven              | Elastic_Net_Cox       |            9 |               0.731 |            0.816 |                 0.773 |                  5 |
| seven              | Elastic_Net_Cox       |            9 |               0.738 |            0.807 |                 0.773 |                  5 |
| seven_interactions | Elastic_Net_Cox       |           14 |               0.737 |            0.806 |                 0.772 |                  5 |
| seven_interactions | Elastic_Net_Cox       |           14 |               0.744 |            0.796 |                 0.770 |                  5 |

## Rescue Headline

Headline rescue model: `Rescue_Cox_PH` selected by TCGA-only composite CV.

## External TNBC Rows vs genefu PAM50 ROR-S

| cohort    | subgroup   | headline_model   | comparator        |   n |   events |   delta_cindex |   ci_low |   ci_high |       p |
|:----------|:-----------|:-----------------|:------------------|----:|---------:|---------------:|---------:|----------:|--------:|
| GSE20685  | tnbc       | Rescue_Cox_PH    | PAM50_RORS_genefu |   0 |        0 |        nan     |  nan     |   nan     | nan     |
| GSE96058  | tnbc       | Rescue_Cox_PH    | PAM50_RORS_genefu | 143 |       26 |          0.199 |    0.078 |     0.328 |   0.001 |
| METABRIC  | tnbc       | Rescue_Cox_PH    | PAM50_RORS_genefu | 320 |      168 |          0.101 |    0.048 |     0.159 |   0.000 |
| TCGA-BRCA | tnbc       | Rescue_Cox_PH    | PAM50_RORS_genefu | 157 |       25 |          0.408 |    0.247 |     0.566 |   0.000 |

## Random-Effects Meta-Analysis

| label                      | cohorts                               | met   |   meta_n_studies |   meta_effect |   meta_se |   meta_ci_low |   meta_ci_high |   meta_z |   meta_p |   meta_tau2 |   meta_i2 |
|:---------------------------|:--------------------------------------|:------|-----------------:|--------------:|----------:|--------------:|---------------:|---------:|---------:|------------:|----------:|
| all_protocol_tnbc_cohorts  | ['GSE96058', 'METABRIC', 'TCGA-BRCA'] | True  |                3 |        0.2225 |    0.0845 |        0.0570 |         0.3881 |   2.6344 |   0.0084 |      0.0178 |    0.8522 |
| external_only_tnbc_cohorts | ['GSE96058', 'METABRIC']              | True  |                2 |        0.1333 |    0.0458 |        0.0435 |         0.2231 |   2.9101 |   0.0036 |      0.0023 |    0.4880 |

## External Validation Metrics for Headline and PAM50

| cohort    | model             | subgroup   |    n |   events |   harrell_c |   harrell_c_ci_low |   harrell_c_ci_high | status         |
|:----------|:------------------|:-----------|-----:|---------:|------------:|-------------------:|--------------------:|:---------------|
| GSE20685  | Rescue_Cox_PH     | overall    |  327 |       83 |       0.707 |              0.650 |               0.760 | ok             |
| GSE20685  | Rescue_Cox_PH     | tnbc       |    0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| GSE20685  | PAM50_RORS_genefu | overall    |  327 |       83 |       0.635 |              0.566 |               0.693 | ok             |
| GSE20685  | PAM50_RORS_genefu | tnbc       |    0 |        0 |     nan     |            nan     |             nan     | too_few_events |
| GSE96058  | Rescue_Cox_PH     | overall    | 3273 |      336 |       0.760 |              0.733 |               0.785 | ok             |
| GSE96058  | Rescue_Cox_PH     | tnbc       |  143 |       26 |       0.776 |              0.694 |               0.854 | ok             |
| GSE96058  | PAM50_RORS_genefu | overall    | 3273 |      336 |       0.646 |              0.615 |               0.674 | ok             |
| GSE96058  | PAM50_RORS_genefu | tnbc       |  143 |       26 |       0.577 |              0.487 |               0.667 | ok             |
| METABRIC  | Rescue_Cox_PH     | overall    | 1979 |     1143 |       0.646 |              0.628 |               0.662 | ok             |
| METABRIC  | Rescue_Cox_PH     | tnbc       |  319 |      168 |       0.626 |              0.586 |               0.672 | ok             |
| METABRIC  | PAM50_RORS_genefu | overall    | 1979 |     1143 |       0.595 |              0.578 |               0.612 | ok             |
| METABRIC  | PAM50_RORS_genefu | tnbc       |  319 |      168 |       0.525 |              0.481 |               0.569 | ok             |
| TCGA-BRCA | Rescue_Cox_PH     | overall    | 1064 |      150 |       0.756 |              0.702 |               0.803 | ok             |
| TCGA-BRCA | Rescue_Cox_PH     | tnbc       |  157 |       25 |       0.848 |              0.752 |               0.922 | ok             |
| TCGA-BRCA | PAM50_RORS_genefu | overall    | 1064 |      150 |       0.615 |              0.560 |               0.665 | ok             |
| TCGA-BRCA | PAM50_RORS_genefu | tnbc       |  157 |       25 |       0.440 |              0.304 |               0.566 | ok             |
