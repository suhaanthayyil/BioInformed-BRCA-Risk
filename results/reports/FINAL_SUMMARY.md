# Final Summary

Generated: 2026-05-16

## Primary Endpoint

Status: **NOT MET**

Locked Phase 6 model: `Gradient_Boosted_Survival`

Official PAM50 comparator: `genefu::rorS` (`PAM50_ROR_official`)

TNBC random-effects delta C-index vs official PAM50-ROR:

- Effect: `+0.0144`
- 95% CI: `[-0.0472, 0.0760]`
- p-value: `0.6466`
- Cohorts: GSE96058, METABRIC, TCGA-BRCA

## Rescue Analysis

Status: **NOT MET**

Selected rescue model: `Rescue_DeepSurv`

Selection policy: TCGA-only composite CV, `0.5 * overall C-index + 0.5 * TNBC C-index`.

Rescue all-protocol TNBC result:

- Effect: `+0.1284`
- 95% CI: `[-0.1312, 0.3880]`
- p-value: `0.3322`

Rescue external-only TNBC result:

- Effect: `+0.0023`
- 95% CI: `[-0.0624, 0.0669]`
- p-value: `0.9452`

## Key Interpretation

The expanded feature/model rescue improved TCGA internal performance but did not produce externally validated TNBC superiority over official PAM50-ROR. The current repository should not claim that the ML pathway model beats PAM50 in TNBC.

## Main Deliverables

- `results/Table_S4_rescue_internal_cv.csv`
- `results/Table_S5_rescue_external_validation.csv`
- `results/Table_S6_rescue_head_to_head.csv`
- `results/Table_S7_rescue_posthoc_top_candidate_screen.csv`
- `results/reports/rescue_analysis_qc.md`
- `data/processed/rescue_analysis.metadata.json`
- `data/processed/rescue_primary_endpoint.json`

## Recommended Submission Framing

Use an honest pivot: competitive and interpretable pathway ML with official clinical baseline comparisons, not TNBC superiority. The strongest remaining story is biological interpretation and transparent negative/competitive benchmarking.
