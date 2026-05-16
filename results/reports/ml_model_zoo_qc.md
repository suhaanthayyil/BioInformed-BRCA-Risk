# Phase 5 ML Model Zoo QC

Generated: 2026-05-16T12:38:29-04:00

Training cohort: TCGA-BRCA, n=213, events=132
Features: Pathway_Immune, Pathway_Proliferation, Pathway_DNA_Repair, Pathway_Metabolism, Pathway_Stromal_EMT, Pathway_Apoptosis_Stress, Pathway_Hormone, age_at_dx, stage_ordinal

## Internal 5-fold CV

| model                     |   mean_cv_cindex |   sd_cv_cindex | fold_cindexes                                                                                        | config                                                                                                          | artifact_path                        |   auc_5y |   auc_10y |   brier_5y |   delta_vs_cox |
|:--------------------------|-----------------:|---------------:|:-----------------------------------------------------------------------------------------------------|:----------------------------------------------------------------------------------------------------------------|:-------------------------------------|---------:|----------:|-----------:|---------------:|
| Gradient_Boosted_Survival |            0.642 |          0.034 | [0.6748554913294798, 0.6324850299401198, 0.6642754662840746, 0.6503215434083601, 0.5871913580246914] | {"n_estimators": 100, "learning_rate": 0.03, "max_depth": 1, "subsample": 1.0}                                  | models/gradient_boosted_survival.pkl |    0.656 |     0.623 |      0.236 |          0.042 |
| Stacked_Ensemble          |            0.626 |        nan     | []                                                                                                   | {"top_models": ["Gradient_Boosted_Survival", "DeepSurv", "Elastic_Net_Cox"]}                                    | models/stacked_ensemble.pkl          |    0.674 |     0.643 |      0.228 |          0.026 |
| DeepSurv                  |            0.609 |          0.039 | [0.638728323699422, 0.5688622754491018, 0.6441893830703013, 0.6286173633440515, 0.5648148148148148]  | {"hidden_dims": [128, 64], "dropout": 0.2, "lr": 0.0005, "weight_decay": 0.0001, "epochs": 150, "patience": 20} | models/deepsurv.pt                   |    0.640 |     0.602 |      0.246 |          0.009 |
| Elastic_Net_Cox           |            0.601 |          0.034 | [0.6069364161849711, 0.5553892215568862, 0.5868005738880918, 0.6495176848874598, 0.6080246913580247] | {"alpha": 0.001, "l1_ratio": 0.1}                                                                               | models/elastic_net_cox.pkl           |    0.653 |     0.556 |      0.244 |          0.001 |
| Random_Survival_Forest    |            0.601 |          0.044 | [0.6517341040462428, 0.5793413173652695, 0.6327116212338594, 0.5980707395498392, 0.5416666666666666] | {"n_estimators": 500, "max_depth": 4, "min_samples_leaf": 1, "max_features": "sqrt"}                            | models/random_survival_forest.pkl    |    0.648 |     0.612 |      0.239 |          0.001 |
| Cox_PH                    |            0.600 |          0.036 | [0.6083815028901735, 0.5553892215568862, 0.5824964131994261, 0.6527331189710611, 0.6003086419753086] | {"penalizer": 0.01}                                                                                             | models/cox_ph.pkl                    |    0.652 |     0.556 |      0.243 |          0.000 |

Headline model by TCGA CV: Gradient_Boosted_Survival (C-index 0.642, delta vs Cox 0.042).

Grid mode: compact local grid. This is explicitly recorded and should not be described as the full protocol grid.
XGBoost Cox is disabled by default because the local native library repeatedly terminated the process without a Python traceback.
