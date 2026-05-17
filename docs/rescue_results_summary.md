# Rescue Analysis Results Summary

Generated: 2026-05-16

## Bottom Line

The transparent rescue run improved TCGA cross-validation performance but did **not** rescue the pre-registered TNBC endpoint against official `genefu` PAM50-ROR.

The locked Phase 6 result remains the primary result:

- Headline locked model: `Gradient_Boosted_Survival`
- TNBC delta C-index vs official PAM50-ROR: `+0.0144`
- 95% CI: `[-0.0472, 0.0760]`
- p-value: `0.6466`
- Endpoint status: `NOT MET`

## Rescue Run

The rescue run evaluated 391 TCGA-selected candidates across expanded pathway feature sets, pathway interactions, TNBC-weighted training variants, RSF, GBSA, Cox/Coxnet, and DeepSurv. XGBoost AFT was disabled by default because local native survival runs terminated without Python traceback.

Selected rescue model:

- Model: `Rescue_DeepSurv`
- Feature set: `component7_interactions`
- Features: 43
- TCGA overall CV C-index: `0.6346`
- TCGA TNBC CV C-index: `0.6475`
- TCGA-only composite CV C-index: `0.6411`

## Rescue Endpoint Result

All protocol TNBC cohorts:

- Cohorts: GSE96058, METABRIC, TCGA-BRCA
- Delta C-index vs official PAM50-ROR: `+0.1284`
- 95% CI: `[-0.1312, 0.3880]`
- p-value: `0.3322`
- Endpoint status: `NOT MET`

External-only TNBC cohorts:

- Cohorts: GSE96058, METABRIC
- Delta C-index vs official PAM50-ROR: `+0.0023`
- 95% CI: `[-0.0624, 0.0669]`
- p-value: `0.9452`
- Endpoint status: `NOT MET`

## Interpretation

The rescue model overfit the TCGA training cohort and did not generalize enough to support a significant-advantage claim. Its TCGA TNBC delta was large (`+0.4244`), but GSE96058 was negative (`-0.0386`) and METABRIC was only slightly positive (`+0.0090`).

The strongest post-hoc external-only screen among near-top TCGA candidates was a `seven_interactions` Random Survival Forest, with mean external-only TNBC delta around `+0.0133`; this is still below the `+0.03` target and was not the locked TCGA-selected rescue headline.

## Submission Implication

Do not submit a claim that pathway ML has a significant TNBC advantage over PAM50-ROR. The defensible path is an honest pivot: pathway-level ML is competitive in some settings, official PAM50 remains difficult to exceed externally, and the stronger contribution should be biological interpretation/subtype-specific pathway findings rather than clinical test advantage.
