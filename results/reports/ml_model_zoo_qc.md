# Phase 5 ML Model Zoo QC

Generated: 2026-10-03T17:59:37-04:00

Training cohort: TCGA-BRCA, n=1064, events=150
Features: Pathway_Immune, Pathway_Proliferation, Pathway_DNA_Repair, Pathway_Metabolism, Pathway_Stromal_EMT, Pathway_Apoptosis_Stress, Pathway_Hormone

## Internal 5-fold CV

| model                     |   mean_cv_cindex |   sd_cv_cindex |   single_split_cv_cindex |   n_cv_folds | fold_cindexes                                                                                        | config                                                                                                          | artifact_path                        |   auc_5y |   auc_10y |   brier_5y |   delta_vs_cox | headline_eligible   |
|:--------------------------|-----------------:|---------------:|-------------------------:|-------------:|:-----------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------|:-------------------------------------|---------:|----------:|-----------:|---------------:|:--------------------|
| Cox_PH                    |            0.614 |          0.056 |                    0.619 |           25 | [0.6524390243902439, 0.6642501132759402, 0.576530612244898, 0.5239423523942353, 0.6776649746192893]  | {"penalizer": 0.1}                                                                                              | models/cox_ph.pkl                    |    0.628 |     0.492 |      0.251 |          0.000 | False               |
| Elastic_Net_Cox           |            0.614 |          0.054 |                    0.619 |           25 | [0.646680216802168, 0.6909832351608518, 0.5663265306122449, 0.5109251510925151, 0.6791153009427121]  | {"alpha": 0.01, "l1_ratio": 0.1}                                                                                | models/elastic_net_cox.pkl           |    0.620 |     0.519 |      0.257 |         -0.000 | False               |
| Stacked_Ensemble          |            0.613 |        nan     |                  nan     |            5 | []                                                                                                   | {"top_models": ["Cox_PH", "Elastic_Net_Cox", "DeepSurv"]}                                                       | models/stacked_ensemble.pkl          |    0.621 |     0.551 |      0.251 |         -0.001 | False               |
| DeepSurv                  |            0.593 |          0.071 |                    0.610 |           25 | [0.6029810298102981, 0.7109198006343452, 0.532574568288854, 0.5178986517898652, 0.6860043509789703]  | {"hidden_dims": [128, 64], "dropout": 0.3, "lr": 0.0003, "weight_decay": 0.0001, "epochs": 300, "patience": 30} | models/deepsurv.pt                   |    0.598 |     0.552 |      0.256 |         -0.021 | True                |
| Random_Survival_Forest    |            0.577 |          0.059 |                    0.589 |           25 | [0.6161924119241192, 0.674218396012687, 0.5302197802197802, 0.5490469549046955, 0.5754169688179841]  | {"n_estimators": 500, "max_depth": 5, "min_samples_leaf": 40, "max_features": 0.5}                              | models/random_survival_forest.pkl    |    0.582 |     0.473 |      0.261 |         -0.037 | True                |
| Gradient_Boosted_Survival |            0.574 |          0.065 |                    0.581 |           25 | [0.5789295392953929, 0.6273221567739012, 0.5592621664050236, 0.5657833565783357, 0.5717911530094271] | {"n_estimators": 100, "learning_rate": 0.01, "max_depth": 1, "subsample": 0.8}                                  | models/gradient_boosted_survival.pkl |    0.547 |     0.485 |      0.276 |         -0.040 | True                |

Headline model (best registered non-linear candidate: RSF, GBSA or DeepSurv) by TCGA CV: DeepSurv (C-index 0.593, delta vs Cox -0.021).

Grid mode: expanded grid (RSF 75, GBSA 81, DeepSurv 40 configurations); top 5 per model re-scored with 5x5 repeated CV; selection uses TCGA only.
XGBoost Cox is disabled by default because the local native library repeatedly terminated the process without a Python traceback.
