# Pre-Registered Primary Endpoint

**Project**: BRCA-PathwayML v2 (revision after Communications Medicine + Scientific Reports desk rejections).

**Target venue**: npj Breast Cancer (primary), BMC Artificial Intelligence (current submitted venue after analysis reframing).

**Primary claim**: A non-linear machine learning model (Random Survival Forest or Gradient-Boosted Survival or DeepSurv, whichever wins on internal CV) trained on 7 interpretable biological pathway scores is non-inferior to PAM50 and the Oncotype DX 21-gene surrogate on overall discrimination AND is superior to both on the triple-negative (TNBC) subgroup, across at least three independent breast cancer cohorts.

**Primary endpoint**: TNBC meta-analyzed delta C-index (best ML model vs PAM50-ROR) across all cohorts with TNBC patients, random-effects meta-analysis.

**Success threshold (PRE-REGISTERED)**: Meta-analyzed delta C-index >= +0.03 in TNBC, p < 0.05.

**Secondary claims**:
(a) ML model outperforms Cox baseline by >= +0.02 C-index on at least one cohort (showing real ML adds value).
(b) SHAP / model-class-appropriate attribution identifies immune-axis pathway as dominant in TNBC and HER2+, recovering known biological pattern.
(c) Combined model (ML pathway + PAM50 + clinical) shows IDI > 0.05 over PAM50 alone.
(d) ML model non-inferior to PAM50 (one-sided, margin -0.02) on overall cohort.

**Pre-registration timestamp**: 2026-05-16T11:10:43-04:00

**Modifications after analysis begins**: NOT PERMITTED.

## Locked-artifact manifest

The weights trained and locked on 2026-05-16 (TCGA-BRCA n = 213, nine features
including a stage covariate that was missing for every external patient) are
kept unchanged in `models/locked_v1/` and fingerprinted in
`docs/MODEL_SHA256.txt`. `scripts/locked_model_sensitivity.py` scores them on
every eligible external patient (Table S19).

The current weights in `models/` come from the revised analysis (all eligible
TCGA-BRCA primary tumours, the seven registered pathway scores only); see `models/README.md` and
`data/processed/ml_model_zoo.metadata.json`.

**Timeline (corrected 2026-10-01).** The registration text above was committed
at `7783f4f` on 2026-05-16 11:11 (-04:00). It was not the first analysis in
this repository. The first commit (`9bef60a`, 2026-04-01) already contained the
earlier two-cohort study that was posted as the April 2026 preprint: TCGA-BRCA
models with a different seven-pathway definition and binary 5-year outcome,
applied to GSE96058/SCAN-B (`results/legacy_v1_preprint/transfer_results.csv`, notebooks 03-14).
`7783f4f` itself added `paper/main.tex`, which reported those GSE96058 results.

What was new after the lock (all on or after 2026-05-16 11:17): METABRIC and
GSE20685 acquisition, cohort harmonization, the MSigDB pathway definitions,
the survival model zoo, genefu PAM50 ROR-S scoring, the TNBC primary endpoint
and every survival result reported in the BMC submission. GSE96058 and TCGA
had been analysed before the lock under the earlier design.

The sentence "Modifications after analysis begins: NOT PERMITTED" was not
followed: this file was edited on 2026-05-24 (venue) and 2026-06-22 (manifest
and timeline note), and the analysis was revised in October 2026 after the
journal's pre-publication checks. Every deviation from the registered plan,
with its reason, is listed in `docs/REGISTRATION_DEVIATIONS.md`.
