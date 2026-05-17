# Story Lock

## Pre-registered claim (LOCKED, NOT MET)

Pathway-based ML exceeds PAM50-ROR in TNBC by >= +0.03 C-index, p < 0.05.

Result: meta delta = +0.0144, 95% CI [-0.0472, 0.0760], p = 0.6466.

Status: NOT MET. This is reported transparently.

## Primary finding

Interpretable pathway-based ML achieves discrimination comparable to official PAM50-ROR in four-cohort breast cancer prognosis (4,532 harmonized samples; 4,003 evaluable in the official PAM50-ROR head-to-head analysis), while providing biological attribution that subtype labels alone cannot. ML provides modest internal gains over Cox (+0.042 C-index) but does not significantly improve upon PAM50-ROR externally.

## Secondary findings (exploratory, post-hoc)

- Within-subtype risk stratification: comparable, with no consistent within-subtype advantage over PAM50-ROR.
- Calibration: headline ML showed lower mean ICI than PAM50-ROR (0.174 vs 0.204).
- Net benefit (DCA): mean 5-year net benefit delta versus PAM50-ROR was +0.0156 across cohorts and thresholds.
- Stability of pathway attributions: mean pairwise Spearman rho=0.395 (95% empirical interval -0.356 to 0.917); top-ranked features: Pathway_Hormone, age_at_dx, Pathway_DNA_Repair.

## What we do NOT claim

- ML has a TNBC advantage over PAM50-ROR. This was pre-registered and not supported.
- The rescue model generalizes. It does not, and it shows a clear overfit pattern.
- A significant advantage over any established clinical genomic test.

## Authors

- Suhaan Thayyil (first author, corresponding)
- Eshaan Nidee (co-author)
