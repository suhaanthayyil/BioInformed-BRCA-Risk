# BRCA-PathwayML

Interpretable pathway-based machine learning for breast cancer prognosis: a four-cohort benchmark against PAM50-ROR.

**Authors:** Suhaan Thayyil, Eshaan Nidee

## Overview

This repository contains the code, analysis pipelines, and results for a four-cohort study comparing a pathway-based survival machine learning model (Gradient Boosted Survival Analysis) against official genefu PAM50-ROR for breast cancer prognosis.

The pre-registered primary endpoint (TNBC meta-analyzed delta C-index >= +0.03, p < 0.05) was **not met**. The observed result was delta = +0.0144, 95% CI [-0.0472, +0.0760], p = 0.6466. All secondary analyses are exploratory.

## Key Results

- **Harmonized database:** 4,532 samples across four cohorts (TCGA-BRCA, GSE96058/SCAN-B, METABRIC, GSE20685).
- **Analysis sets:** 4,003 patients in patient-characteristics summaries; 4,002 metric-evaluable patients in the official PAM50-ROR discrimination analysis.
- **Headline model:** Gradient Boosted Survival (internal TCGA CV C-index = 0.642, delta vs Cox = +0.042).
- **TNBC primary endpoint:** NOT MET (meta delta = +0.0144, p = 0.6466).
- **Overall discrimination:** comparable to PAM50-ROR across external cohorts.

## Submitted Manuscript Crosswalk

The submitted BMC Artificial Intelligence manuscript is mapped to repository files in `docs/SUBMITTED_MANUSCRIPT_REPOSITORY_CROSSWALK.md`.

Reviewer-facing traceability helpers:

- `scripts/check_submission_consistency.py` verifies the central submitted values against repository CSV/JSON outputs.
- `results/Table_S2_internal_cv_manuscript_crosswalk.csv` maps internal cross-validation values between the manuscript summary and exact reproducible output.
- `results/Table_1_denominator_crosswalk.csv` documents patient-characteristics versus metric-evaluable denominators.

## Repository Structure

```
src/              Python modules (models, survival analysis, data loading, features)
src/ml/           DeepSurv implementation
R/                Official genefu PAM50-ROR scoring (pam50.R)
scripts/          Analysis pipeline scripts
notebooks/        Jupyter notebooks (01-14: original analyses)
results/          All results tables and QC reports
figures/          Generated figures
paper/            Reviewer-facing LaTeX manuscript source and DOI-backed bibliography
data/             Raw, processed, and clinical data
models/           Trained model artifacts (gitignored)
tests/            Pytest test suite
docs/             Pre-registration endpoint, story lock, cohort notes
```

## Manuscript LaTeX

The repository includes a reviewer-facing LaTeX manuscript source at `paper/main.tex`, a DOI-backed BibTeX bibliography at `paper/main.bib`, a compiled PDF at `paper/main.pdf`, and a Word export at `paper/main.docx`. The LaTeX source mirrors the submitted competitive-not-superior framing and is intended for traceability, not as a claim that the manuscript source is the journal's copyedited version.

Build commands:

```bash
cd paper && tectonic --keep-logs main.tex
cd .. && python3 scripts/build_manuscript_docx.py
```

## Data Sources

- **TCGA-BRCA:** GDC Data Portal (https://portal.gdc.cancer.gov)
- **GSE96058/SCAN-B:** NCBI GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058)
- **METABRIC:** cBioPortal (https://www.cbioportal.org/study/summary?id=brca_metabric)
- **GSE20685:** NCBI GEO (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20685)

Raw data are not redistributed in this repository.

## Pre-registration

The primary endpoint was locked before any external validation was run. See `docs/PRIMARY_ENDPOINT.md` (committed at `7783f4f0bfaaa6bdc611c78d33ccda621c6b243d`, before any results commits).

## Reproduction

1. Install Python dependencies: `pip install -r requirements.txt`
2. Install R with the `genefu` package for official PAM50-ROR scoring.
3. Pipeline scripts are in `scripts/`. Key steps:
   - `scripts/harmonize_cohorts.py` -- build unified cohort database
   - `scripts/score_pathways.py` -- compute pathway features
   - `scripts/train_ml_zoo.py` -- train survival models
   - `scripts/external_validation.py` -- run external validation
   - `scripts/compute_clinical_baselines.py` -- compute clinical genomic baselines
   - Phase 3 exploratory scripts: `within_subtype_analysis.py`, `calibration_analysis.py`, `dca_analysis.py`, `stability_analysis.py`
   - Sensitivity analyses: `external_only_tnbc_sensitivity.py`, `feature_ablation.py`, `pathway_scoring_sensitivity.py`

## License

MIT License. See `LICENSE`.
