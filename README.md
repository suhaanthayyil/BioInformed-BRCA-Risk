# BRCA-PathwayML

BMC Cancer submission package for an interpretable pathway-based machine learning benchmark in breast cancer prognosis.

**Authors:** Suhaan Thayyil; Eshaan Nidee

## Current Framing

This repository contains a four-cohort analysis comparing a locked pathway-based survival machine learning model against official genefu PAM50-ROR. The pre-registered TNBC endpoint was not met. The active BMC Cancer framing is transparent benchmarking with comparable discrimination, modest secondary calibration and decision-curve signals, and clear limitations.

## Key Locked Result

TNBC delta C-index for Gradient Boosted Survival versus PAM50-ROR = +0.0144, 95% CI [-0.0472, +0.0760], p=0.6466. The pre-registered threshold was delta >= +0.03 with p < 0.05. It was not met.

## Reproducibility

- Analysis cohort: 4003 patients with survival and official PAM50-ROR across four cohorts. Harmonized database contains 4532 samples.
- Main results: `results/`
- BMC submission files: `deliverables/`
- Locked endpoint: `docs/PRIMARY_ENDPOINT.md`
- Phase Three story lock: `docs/STORY.md`

## Data Sources

TCGA-BRCA, GSE96058/SCAN-B, METABRIC from cBioPortal, and GSE20685. Raw data are not redistributed here.

## License

MIT License. See `LICENSE`.
