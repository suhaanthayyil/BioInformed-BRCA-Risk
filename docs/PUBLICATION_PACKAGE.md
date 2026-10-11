# Publication Package Evidence Map

Release `v1.2.1` is the code and reproducibility package for the revised 6,644-patient analysis. Manuscript, Supplementary Information, and journal correspondence files are kept in the private journal-submission package and are not distributed through GitHub.

## Registration timing

The immutable registration text is reproduced in `docs/PRIMARY_ENDPOINT.md`. Its corrected timeline distinguishes the April 2026 TCGA-BRCA-to-GSE96058 binary 5-year analysis from the analyses first performed after registration: acquisition and use of METABRIC and GSE20685, the current MSigDB pathway definitions, the time-to-event survival model zoo, genefu PAM50 ROR-S scoring, the TNBC meta-analysis, and the other survival analyses reported in the BMC manuscript.

Evidence:

- `docs/PRIMARY_ENDPOINT.md`
- `docs/REGISTRATION_DEVIATIONS.md`
- `results/Table_S20_registration_deviations.csv`
- Supplementary Section S3 and Table S20

## Registered and reported analyses

Supplementary Table S18 contains all eight registered endpoints. This includes both Oncotype DX surrogate dispositions: the overall ER-positive comparison and the TNBC comparison reported as not evaluable. Supplementary Table S20 contains all 16 registered-versus-reported changes and their reasons. Tables S18-S21 are embedded in the journal-submission Supplementary Information; their public machine-readable sources are:

- `results/Table_S18_registered_endpoints.csv`
- `results/Table_S19_locked_model_full_cohorts.csv`
- `results/Table_S20_registration_deviations.csv`
- `results/Table_S21_age_adjusted_comparison.csv`

## Cohort selection and event counts

`results/Figure1_cohort_flow_counts.json` records source counts, exclusions, analysed patients, events, and PAM50-evaluable patients. `results/Table_1_denominator_crosswalk.csv` independently cross-checks the table and metric denominators. The totals are 6,644 analysed and 6,643 evaluable against PAM50 ROR-S.

## Clinical variables and stage

The primary model consumes the seven entries in `src.clinical.PRIMARY_FEATURES`; age and stage are excluded. `src/clinical.stage_to_ordinal` parses TCGA pathologic stage, numeric METABRIC `TUMOR_STAGE`, and the GSE20685 TNM strings created during ingestion. GSE96058 is assigned no stage because its public release does not provide it. Stage-containing sensitivity analyses in `scripts/feature_ablation.py` train and evaluate only rows with recorded stage.

The committed unified cohort database contains parsed stage for 1,042 TCGA-BRCA, 1,764 METABRIC, and 325 GSE20685 records, and zero GSE96058 records. Tests cover representative formats from every contributing cohort in `tests/test_clinical.py`.

## Verification

Run the repository-only publication gate:

```bash
python scripts/check_publication_readiness.py
```

When working from the private submission workspace, include the Supplementary Information in the gate:

```bash
python scripts/check_publication_readiness.py \
  --supplement _submission/BMC_Supplement_prepublication_revision_v1.2.1.docx
```

The broader reproducibility gates are:

```bash
make test
make consistency
```
