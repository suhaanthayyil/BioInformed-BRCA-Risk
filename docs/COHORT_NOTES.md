# Cohort Notes

Last updated: 2026-05-16T11:16:00-04:00

## TCGA-BRCA

- Existing repository files: `data/01_tcga_expression_normalized.csv`, `data/clinical/01_tcga_clinical.csv`, and `data/processed/02_tcga_feature_matrix.csv`.
- Current usable clinical rows: 213 patients with `time_to_event`, `event_status`, `high_risk`, and published `pam50_subtype` fields.
- Expression normalization in the v1 repo is already preprocessed/log-normalized; v2 harmonization must not assume that TCGA centering/scaling transfers to external cohorts.

## SCAN-B / GSE96058

- Source: NCBI GEO supplemental file `GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz`.
- Download SHA256: `3c717baf7960e1f1477a72744a399a42fabba113c9e7ecb3d55997574d7e9732`.
- Parsed output: `data/processed/02_gse96058_expression.parquet`.
- Parsed shape: 30,865 gene rows by 3,409 sample columns, with expression in gene-row/sample-column orientation.
- Clinical rows currently available in-repo: 1,483 rows in `data/clinical/01_gse96058_clinical.csv`.
- Normalization note: the downloaded expression file is already transformed by the SCAN-B/GEO release. Cross-cohort analyses should compute pathway scores within cohort or use rank/standardization procedures that do not leak training cohort scaling into validation cohorts.

## METABRIC

- Primary protocol URL `https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz` returned HTTP 403 on 2026-05-16.
- Current official cBioPortal asset endpoint used: `https://datahub.assets.cbioportal.org/brca_metabric.tar.gz`.
- Archive SHA256: `6d4683477d6b37a2d7edbedc0df610f67bc456f99e5e1bef6219f37b633a55f7`.
- Extracted cBioPortal files include `data_clinical_patient.txt`, `data_clinical_sample.txt`, `data_mrna_illumina_microarray.txt`, `data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt`, and `data_mutations.txt`.
- Current DataHub filenames differ from older protocol names: `data_mrna_illumina_microarray_zscores_ref_diploid_samples.txt` replaces the expected older `data_mRNA_median_Zscores.txt`, and `data_mutations.txt` replaces `data_mutations_extended.txt`.
- Parsed outputs: `data/processed/03_metabric_expression.parquet` and `data/processed/03_metabric_clinical.parquet`.
- Parsed expression shape: 20,603 gene rows by 1,980 sample columns, cBioPortal z-scores relative to diploid samples.
- Parsed clinical shape: 2,509 sample rows, with a first-pass TNBC flag count of 320.
- Normalization note: METABRIC expression z-scores are not directly comparable to RNA-seq transformed expression values. V2 models must treat METABRIC as an external validation cohort and avoid training-time normalization leakage.

## Pending Data

- MSigDB gene sets are intentionally not downloaded until explicit license approval at STOP GATE 1.
- Optional fourth cohort acquisition (GSE20685 or similar) has not started because the protocol stops at the MSigDB license gate first.

