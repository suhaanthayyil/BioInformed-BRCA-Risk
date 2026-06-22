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

The trained weights that produced the locked external-validation results are
committed under `models/` and fingerprinted in `docs/MODEL_SHA256.txt`
(`shasum -a 256 -c docs/MODEL_SHA256.txt`). The headline model is
`gradient_boosted_survival.pkl`; its features, training cohort (TCGA-BRCA,
n = 213), seed (42), and selected configuration are recorded in
`models/README.md` and `data/processed/ml_model_zoo.metadata.json`. This lets a
reviewer confirm the scored weights were not altered after the endpoint lock.

**Timeline note.** This pre-registration text was committed at
`7783f4f0bfaaa6bdc611c78d33ccda621c6b243d`, before any results commits. The
internal pre-registration timestamp above (2026-05-16) precedes the first
public Git commit date (2026-05-17) because the endpoint was fixed locally
before the repository history was published; the ordering of the lock relative
to results is preserved in the Git history.
