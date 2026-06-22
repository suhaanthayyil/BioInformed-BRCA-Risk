# Data Dictionary (Codebook)

This codebook defines the variables used in the BRCA-PathwayML four-cohort
benchmark. Every column listed below was read directly from the committed
repository files. Columns are tagged by role:

- **ID** — identifier / join key (not a model variable)
- **INPUT** — model feature consumed by the survival ML models
- **OUTCOME** — survival endpoint (time and/or event)
- **BASELINE** — comparator output (PAM50-ROR), not a model input
- **AUX** — auxiliary / descriptive field (subtype, receptor status, treatment, etc.)

The headline survival model (Gradient Boosted Survival) is trained on a
fixed 9-column feature schema: the seven harmonized pathway features plus two
clinical covariates (`age_at_dx`, `stage_ordinal`). Source: `scripts/train_ml_zoo.py`
(`PATHWAY_FEATURES`, `feature_columns`) and `data/processed/ml_model_zoo.metadata.json`
(`features`).

---

## 1. The seven harmonized pathway features

Defined in `src/pathways.py` (`SEVEN_PATHWAY_COMPONENTS`) and materialized in
`data/processed/04_pathway_features.parquet`. Each feature is the **mean of the
z-scored, rank-based ssGSEA-style scores** of its component MSigDB v2024.1.Hs
gene sets, computed **within cohort** (rank percentile within sample, then
gene-set scores z-scored across samples within cohort, then averaged across the
component sets). All seven are continuous, dimensionless, and approximately
zero-centered within each cohort. All are **model INPUTS**.

| Feature (column) | Role | Aggregates these MSigDB gene sets (per `src/pathways.py`) |
|---|---|---|
| `Pathway_Immune` | INPUT | HALLMARK_INTERFERON_GAMMA_RESPONSE, HALLMARK_INTERFERON_ALPHA_RESPONSE, HALLMARK_INFLAMMATORY_RESPONSE, HALLMARK_ALLOGRAFT_REJECTION, HALLMARK_COMPLEMENT |
| `Pathway_Proliferation` | INPUT | HALLMARK_E2F_TARGETS, HALLMARK_G2M_CHECKPOINT, HALLMARK_MYC_TARGETS_V1, HALLMARK_MYC_TARGETS_V2, HALLMARK_MITOTIC_SPINDLE |
| `Pathway_DNA_Repair` | INPUT | HALLMARK_DNA_REPAIR, REACTOME_HDR_THROUGH_HOMOLOGOUS_RECOMBINATION_HRR, REACTOME_NUCLEOTIDE_EXCISION_REPAIR, REACTOME_MISMATCH_REPAIR, KEGG_MEDICUS_REFERENCE_HOMOLOGOUS_RECOMBINATION, KEGG_MEDICUS_REFERENCE_MISMATCH_REPAIR |
| `Pathway_Metabolism` | INPUT | HALLMARK_OXIDATIVE_PHOSPHORYLATION, HALLMARK_GLYCOLYSIS, HALLMARK_FATTY_ACID_METABOLISM, KEGG_MEDICUS_REFERENCE_GLYCOLYSIS |
| `Pathway_Stromal_EMT` | INPUT | HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION, HALLMARK_ANGIOGENESIS, HALLMARK_HEDGEHOG_SIGNALING |
| `Pathway_Apoptosis_Stress` | INPUT | HALLMARK_APOPTOSIS, HALLMARK_P53_PATHWAY, HALLMARK_HYPOXIA |
| `Pathway_Hormone` | INPUT | HALLMARK_ESTROGEN_RESPONSE_EARLY, HALLMARK_ESTROGEN_RESPONSE_LATE, HALLMARK_ANDROGEN_RESPONSE |

Scoring provenance (gene-set counts, coverage, cohort shapes) is in
`data/processed/pathway_scoring.metadata.json` and
`results/reports/pathway_scoring_coverage.csv`. A pathway component is only
included for a cohort if at least 5 of its genes are present (`min_genes=5` in
`rank_ssgsea_scores`); per-cohort component presence is recorded by
`aggregate_seven_pathways`.

---

## 2. Clinical covariates (model inputs)

| Variable | Role | Definition / encoding |
|---|---|---|
| `age_at_dx` | INPUT | Age at diagnosis in years (continuous). Coerced to numeric in `scripts/train_ml_zoo.py`. Sourced from the harmonized `samples` table (`age_at_dx`, DOUBLE). |
| `stage_ordinal` | INPUT | Ordinal AJCC/TNM stage derived from the harmonized `stage_tnm` string by `stage_to_ordinal()` in `scripts/train_ml_zoo.py`: contains "IV" -> 4.0; else "III" -> 3.0; else "II" -> 2.0; else "I" -> 1.0; otherwise NaN. Missing values are median-imputed inside the model pipeline (`SimpleImputer(strategy="median")`). |

---

## 3. Survival outcomes (harmonized)

The modeling pipeline reads survival from the harmonized `survival` table in
`data/processed/unified_cohorts.duckdb` (joined on `sample_id`). These are the
authoritative outcome variables used for all C-index/calibration/DCA analyses.

| Variable | Role | Units / encoding |
|---|---|---|
| `os_days` | OUTCOME (time) | Overall-survival follow-up time in **days** (continuous; DOUBLE). Verified ranges by cohort: TCGA-BRCA ~159-7229, GSE96058 ~56-2474, METABRIC 0-10811, GSE20685 ~146-5150. Training excludes rows with `os_days <= 0` or missing (`scripts/train_ml_zoo.py`). |
| `os_event` | OUTCOME (event) | Overall-survival event indicator: **0 = censored, 1 = death** (BIGINT). |
| `dfs_days` | OUTCOME (time) | Disease-/recurrence-free survival time in **days** (DOUBLE). Populated only for METABRIC (n=2388) and a small subset of GSE20685 (n=25); NULL for TCGA-BRCA and GSE96058. Used by the METABRIC recurrence sensitivity analysis. |
| `dfs_event` | OUTCOME (event) | Disease-/recurrence-free event indicator: **0 = censored, 1 = recurrence/relapse** (BIGINT). |
| `last_followup_days` | AUX | Last-contact time in days (DOUBLE), where available. |

Note on per-cohort raw clinical CSVs: the un-harmonized files
`data/clinical/01_tcga_clinical.csv` and `01_gse96058_clinical.csv` carry their
own `time_to_event` column expressed in **months** (e.g., TCGA 32.2, GSE96058
78.9). Do **not** mix these with the harmonized `os_days` (days). Modeling uses
the harmonized day-scale `survival.os_days` only.

---

## 4. PAM50-ROR comparator (baseline outputs)

Materialized in `data/processed/baselines_pam50.parquet` (long format, one row
per `cohort` x `sample_id` x `baseline`). Two baselines are present:
`PAM50_ROR_official` (genefu, `method_label = genefu_rorS`) and
`PAM50_ROR_surrogate` (`published_subtype_plus_proliferation`). The
**official** baseline is the head-to-head comparator in the manuscript.

| Column | Role | Definition |
|---|---|---|
| `cohort` | ID | Cohort label (TCGA-BRCA, GSE96058, METABRIC, GSE20685). |
| `sample_id` | ID | Harmonized sample identifier (join key). |
| `baseline` | ID | `PAM50_ROR_official` or `PAM50_ROR_surrogate`. |
| `score` | BASELINE | Raw ROR/risk score on the method's native scale. |
| `score_0_100` | BASELINE | ROR score rescaled to **0-100** (verified min 0.0, max 100.0). This is the discrimination comparator score. |
| `method_label` | AUX | `genefu_rorS` (official) or `published_subtype_plus_proliferation` (surrogate). |
| `status` | AUX | `ok` (score computed) or `unavailable` (insufficient inputs). |
| `note` | AUX | Free-text provenance / unavailability reason. |
| `pam50_subtype_genefu` | BASELINE/AUX | genefu-called intrinsic subtype: one of {Basal, Her2, LumA, LumB, Normal}. |
| `rorS_risk` | BASELINE/AUX | genefu ROR-S risk group: {Low, Intermediate, High}. |
| `n_pam50_genes_input` | AUX | Number of PAM50 genes available to the scorer for that sample. |
| `er_status` | AUX | ER status used for risk-group cutpoints. |
| `tnbc_flag` | AUX | Triple-negative indicator (used to define the pre-registered TNBC subgroup). |
| `intrinsic_subtype_pam50_published` | AUX | Published/source-provided intrinsic subtype, for cross-check. |

---

## 5. Per-file column tables

### 5a. `data/clinical/01_tcga_clinical.csv` (TCGA-BRCA, raw)

This file is the published TCGA-BRCA clinical export (UCSC Xena / nature2012
fields) plus six appended harmonized columns. Only the appended/relevant
columns are tabulated; the file additionally contains ~190 upstream TCGA
clinical/`_GENOMIC_ID_*` fields that are not used by the modeling pipeline.

| Column | Role | Notes |
|---|---|---|
| `sampleID` / `sample_id` | ID | TCGA sample barcode / harmonized id. |
| `time_to_event` | OUTCOME (time) | **Months** (raw cohort scale; not the harmonized day-scale used for modeling). |
| `event_status` | OUTCOME (event) | 0 = censored, 1 = death (counts: 0->81, 1->132). |
| `high_risk` | AUX/label | Binary high-risk label present in the raw file (0->139, 1->74). Derived upstream; the modeling pipeline does not consume it as a feature. |
| `stage` | AUX | Source stage string (e.g., "Stage IIA"); harmonized to `stage_ordinal` for modeling. |
| `pam50_subtype` / `PAM50Call_RNAseq` / `PAM50_mRNA_nature2012` | AUX | Published PAM50 subtype calls (LumA, LumB, Basal, Her2, Normal). |
| `OS_Time_nature2012`, `OS_event_nature2012`, `Vital_Status_nature2012` | AUX | Upstream survival fields (superseded by the harmonized `survival` table). |
| `ER_Status_nature2012`, `PR_Status_nature2012`, `HER2_Final_Status_nature2012` | AUX | Receptor status. |
| `Age_at_Initial_Pathologic_Diagnosis_nature2012` / `age_at_initial_pathologic_diagnosis` | AUX | Source age (harmonized to `age_at_dx`). |

### 5b. `data/clinical/01_gse96058_clinical.csv` (GSE96058 / SCAN-B, raw)

| Column | Role | Notes |
|---|---|---|
| `sample_id` | ID | GEO sample id. |
| `scan-b_external_id` | ID | SCAN-B external id. |
| `age_at_diagnosis` | AUX | Source age (harmonized to `age_at_dx`). |
| `overall_survival_days` | OUTCOME (time) | OS time in **days** in this raw file (e.g., 2367, 2168). |
| `overall_survival_event` | OUTCOME (event) | 0 = censored, 1 = death (0->1161, 1->322). |
| `time_to_event` | OUTCOME (time) | **Months** (raw cohort scale; not used for modeling). |
| `event_status` | OUTCOME (event) | 0 = censored, 1 = death (matches `overall_survival_event`). |
| `high_risk` | AUX/label | Binary high-risk label (0->1188, 1->295). Not a model feature. |
| `pam50_subtype` | AUX | Published PAM50 subtype. |
| `er_status`, `pgr_status`, `her2_status`, `ki67_status`, `nhg` | AUX | Clinicopathology / receptor status. |
| `*_prediction_mgc`, `*_prediction_sgc` | AUX | SCAN-B multi-/single-gene classifier predictions. |
| `tumor_size`, `lymph_node_group`, `lymph_node_status` | AUX | Tumor burden / nodal status. |
| `endocrine_treated`, `chemo_treated` | AUX | Treatment flags. |

### 5c. `data/processed/04_pathway_features.parquet` (all cohorts; shape 5929 x 9)

| Column | Role | Notes |
|---|---|---|
| `sample_id` | ID | Harmonized sample id (join key). |
| `cohort` | ID | Cohort label. |
| `Pathway_Immune` | INPUT | See Section 1. |
| `Pathway_Proliferation` | INPUT | See Section 1. |
| `Pathway_DNA_Repair` | INPUT | See Section 1. |
| `Pathway_Metabolism` | INPUT | See Section 1. |
| `Pathway_Stromal_EMT` | INPUT | See Section 1. |
| `Pathway_Apoptosis_Stress` | INPUT | See Section 1. |
| `Pathway_Hormone` | INPUT | See Section 1. |

(`age_at_dx` and `stage_ordinal` are joined in from the harmonized `samples`
table at modeling time, not stored in this parquet.)

### 5d. `data/processed/baselines_pam50.parquet` (all cohorts; shape 9064 x 14)

See Section 4 for the full column table. 4,532 rows per baseline
(`PAM50_ROR_official` and `PAM50_ROR_surrogate`).

---

## 6. Harmonized `samples` and `survival` tables (`unified_cohorts.duckdb`)

The modeling frame in `scripts/train_ml_zoo.py` is built by joining
`04_pathway_features.parquet` to these two tables on `sample_id` / `cohort`.

**`samples`** (per sample): `sample_id` (ID), `cohort` (ID), `patient_id` (ID),
`age_at_dx` (INPUT), `sex`, `race`, `ethnicity`, `stage_tnm` (source of
`stage_ordinal`), `grade`, `er_status`, `pr_status`, `her2_status`,
`intrinsic_subtype_pam50_published`, `tnbc_flag` (BOOLEAN; defines TNBC
subgroup), `primary_site`, `treatment_summary` (all AUX unless noted).

**`survival`** (per sample): `sample_id` (ID), `os_days` (OUTCOME), `os_event`
(OUTCOME), `dfs_days` (OUTCOME), `dfs_event` (OUTCOME), `last_followup_days`
(AUX). See Section 3 for encodings.

---

## 7. Cohort coverage (from `pathway_scoring.metadata.json`)

| Cohort | Samples with pathway features | Role |
|---|---|---|
| TCGA-BRCA | 213 | Training cohort (internal 5-fold CV; headline model fit here) |
| GSE96058 / SCAN-B | 3409 | External validation |
| METABRIC | 1980 | External validation (also DFS/recurrence endpoint) |
| GSE20685 | 327 | External validation (OS only; no receptor status, cannot enter TNBC endpoint) |

Sample counts entering the official PAM50-ROR head-to-head and patient-
characteristics analyses differ from these scoring counts because of metric-
evaluability filters (non-positive survival time, model-specific missingness);
see `results/Table_1_denominator_crosswalk.csv` and
`docs/SUBMITTED_MANUSCRIPT_REPOSITORY_CROSSWALK.md`.
