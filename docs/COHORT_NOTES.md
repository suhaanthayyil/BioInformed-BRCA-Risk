# Cohort Notes

Last updated: 2026-10-03 (cohort rebuild); original notes 2026-05-16.

## Inclusion criteria (revised 2026-10-01)

Every cohort keeps every eligible patient regardless of follow-up length;
censored patients enter the survival analysis at their last follow-up. Counts
per step are written by `scripts/ingest_raw.py` to
`data/processed/ingest_flow.json` and drawn in Figure 1
(`scripts/build_figure1_consort.py`).

## TCGA-BRCA

- Source: UCSC Xena `TCGA.BRCA.sampleMap/HiSeqV2` (log2(RSEM normalized count + 1)), `BRCA_clinicalMatrix`, and TCGA-CDR survival (`survival/BRCA_survival.txt`).
- Inclusion: primary tumour samples (barcode type 01), one per patient, female, not redacted in TCGA-CDR, OS status and OS time > 0 available. 1,097 primary-tumour patients; 14 male or sex missing, 5 redacted and 14 without usable OS excluded; 1,064 analysed, 150 deaths.
- Stage: `pathologic_stage` (AJCC), mapped to I-IV; "Stage X" and discrepant entries treated as missing.
- Receptor status: ER and PR from IHC fields; HER2 positive if IHC or ISH positive, negative if ISH negative or IHC negative.
- Legacy note: the submitted analysis used a 213-sample set inherited from the v1 binary-outcome study (`data/legacy_v1_preprint/01_tcga_clinical.csv`). It kept only patients who died or had at least 5 years of follow-up, and included 37 solid-tissue-normal (-11) and 3 metastatic (-06) samples alongside the primary tumours of the same patients. Its sample IDs are listed in `data/clinical/legacy_v1_cohort_ids.csv`.

## SCAN-B / GSE96058

- Source: GEO supplementary expression file (3,273 patients + 136 technical replicates) and the two series matrices (GPL11154, GPL18573) for clinical data.
- Inclusion: all non-replicate samples with OS; 3,273 analysed, 336 deaths.
- Stage is not released in GEO; tumour size and nodal group are, but no M category, so stage is treated as unavailable.
- Legacy note: the submitted analysis used 1,483 samples from the v1 study (`data/legacy_v1_preprint/01_gse96058_clinical.csv`). That set excluded patients censored before 5 years and contained 60 technical replicates.

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

- Stage: `TUMOR_STAGE` (1-4) used directly; stage 0 (in situ) and missing values treated as missing. 2,509 clinical records, 1,980 with expression; the 529 without an expression profile cannot be scored.

## Pending Data

- MSigDB v2024.1.Hs approval was received on 2026-05-16. The Broad public release directory is `2024.1.Hs/` while filenames carry `v2024.1.Hs`.
- Downloaded and validated GMT files: Hallmark, Reactome, KEGG MEDICUS, PID, and GO BP. Provenance is recorded in `data/processed/msigdb_v2024_1_Hs.metadata.json`.

## GSE20685

- Source: GEO series matrix `GSE20685_series_matrix.txt.gz`, platform GPL570.
- Series matrix SHA256: `e818a4d5834e20bbcf01de515caf856c394abb36a02b9cccc95be05c22d0a279`.
- GPL570 annotation SHA256: `d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394`.
- Parsed outputs: `data/processed/04_gse20685_expression.parquet` and `data/processed/04_gse20685_clinical.parquet`.
- Parsed expression shape: 22,189 gene-symbol rows by 327 sample columns after GPL570 probe-to-symbol mapping and duplicate-symbol averaging.
- Parsed clinical shape: 327 rows with OS event and follow-up duration from GEO sample characteristics.
- Limitation: ER/PR/HER2 receptor status was not available in the parsed series metadata, so GSE20685 can support overall survival validation but cannot contribute to the TNBC primary endpoint unless additional metadata is found.
- Stage: AJCC anatomic stage I-IV derived from the GEO T, N and M fields (`src/clinical.py::tnm_to_stage`).
