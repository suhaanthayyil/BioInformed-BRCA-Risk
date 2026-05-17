# Results Bullets

## Patient cohorts
- Harmonized database: n = 4532 across 4 cohorts.
- Final official PAM50-ROR head-to-head analysis: n = 4003.
- Overall survival events in final analysis set: 1681.
- See Figure 1 and Table 1.

## Discrimination
- Headline model: Gradient Boosted Survival.
- Overall meta delta C-index versus PAM50-ROR: +0.0193, 95% CI [-0.0576, +0.0962], p=0.6223.
- Meta C-index estimate: pathway ML 0.614; PAM50-ROR 0.595; Cox baseline 0.595.
- Interpretation: comparable, not statistically different from PAM50-ROR.
- See Table 2, Table 3, Figure 2 Panel A.

## TNBC subgroup, locked primary endpoint
- TNBC evaluable cohorts: 3.
- TNBC total n in evaluable cohorts: 408; events: 214.
- Locked delta C-index versus PAM50-ROR: +0.0144, 95% CI [-0.0472, +0.0760], p=0.6466.
- Pre-registered success threshold: delta >= +0.03 and p < 0.05.
- Result: NOT MET.
- See Figure 2 Panel B.

## Within-subtype analysis, exploratory
- Basal: delta +0.0349, 95% CI [-0.0123, +0.0822], p=0.147
- Her2: delta -0.0118, 95% CI [-0.1422, +0.1186], p=0.859
- LumA: delta +0.0121, 95% CI [-0.0545, +0.0787], p=0.722
- LumB: delta +0.0423, 95% CI [-0.0742, +0.1588], p=0.477
- Normal: delta -0.0760, 95% CI [-0.2023, +0.0503], p=0.238
- Verdict: no subtype had a consistent meta-analytic advantage meeting the Phase Three reporting rule.
- See Figure 3 and Tables S8-S9.

## Calibration, exploratory
- Mean ICI across cohorts and horizons: pathway ML 0.174; PAM50-ROR 0.204.
- Lower ICI favors the pathway ML model in this exploratory summary.
- See Figure 4 and Table S10.

## Decision curve analysis, exploratory
- Mean 5-year net benefit delta versus PAM50-ROR across cohorts and thresholds: +0.0156.
- Mean categorical NRI at 0.20 threshold: +0.1555.
- Mean IDI across 0.10 and 0.20 thresholds: +0.0707.
- See Figure 5 and Tables S11-S12.

## Stability, exploratory
- Mean pairwise Spearman rho across 100 subsample feature-rank lists: 0.395, empirical interval [-0.356, 0.917].
- Top mean-rank features: Pathway_Hormone, age_at_dx, Pathway_DNA_Repair.
- Limitation: feature-rank stability was below 0.6 and should be discussed as weak.
- See Table S13 and Additional file 5.

## ML versus Cox
- Internal TCGA CV C-index for Gradient Boosted Survival: 0.642.
- Internal delta versus Cox: +0.0420.
- External claim should remain modest because the overall head-to-head delta versus PAM50-ROR was not significant.
