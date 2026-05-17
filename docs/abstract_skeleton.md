# Structured Abstract Skeleton

## Background
[TODO Suhaan: write Background.]

## Methods
- Four cohorts: TCGA-BRCA, GSE96058/SCAN-B, METABRIC, GSE20685.
- Final head-to-head analysis set: n = 4003 with survival and official genefu PAM50-ROR.
- Model: locked Gradient Boosted Survival trained on TCGA-BRCA and externally evaluated without refitting.
- Comparator: official genefu PAM50-ROR.

## Results
- Overall delta C-index versus PAM50-ROR: +0.0193, 95% CI [-0.0576, +0.0962], p=0.6223.
- TNBC primary endpoint: +0.0144, 95% CI [-0.0472, +0.0760], p=0.6466. Not met.
- Calibration: mean ICI 0.174 for pathway ML and 0.204 for PAM50-ROR.
- DCA: mean net benefit delta +0.0156 across cohorts and thresholds.
- Stability: mean pairwise Spearman rho 0.395.

## Conclusions
[TODO Suhaan: write Conclusions with comparable performance and negative primary endpoint.]
