# Data Dictionary (Codebook)

This codebook defines the variables used in the BRCA-PathwayML four-cohort
benchmark. Every column listed below was read directly from the committed
repository files. Columns are tagged by role:

- **ID** — identifier / join key (not a model variable)
- **INPUT** — model feature consumed by the survival ML models
- **OUTCOME** — survival endpoint (time and/or event)
- **BASELINE** — comparator output (PAM50-ROR), not a model input
- **AUX** — auxiliary / descriptive field (subtype, receptor status, treatment, etc.)

The survival models are trained on the seven harmonized pathway features only,
as registered. Source: `src/clinical.py` (`PRIMARY_FEATURES`) and
`data/processed/ml_model_zoo.metadata.json` (`features`). `age_at_dx` and
`stage_ordinal` enter only the sensitivity feature sets
(`scripts/feature_ablation.py`) and the registered combined model (secondary
endpoint c). GSE96058 does not release stage.

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

## 2. Clinical covariates (sensitivity analyses and combined model)

| Variable | Role | Definition / encoding |
|---|---|---|
| `age_at_dx` | SENSITIVITY INPUT | Age at diagnosis in years (continuous). Coerced to numeric in `scripts/train_ml_zoo.py`. Sourced from the harmonized `samples` table (`age_at_dx`, DOUBLE). |
| `stage_ordinal` | SENSITIVITY INPUT | AJCC anatomic stage I-IV from `src/clinical.py::stage_to_ordinal`: Roman-numeral stage strings (TCGA `pathologic_stage`, sub-stages A-C collapsed), numeric 1-4 (METABRIC `TUMOR_STAGE`; 0 = in situ treated as missing) and TNM codes such as `T2.0N1M0` (GSE20685, converted by `tnm_to_stage`). "Stage X", discrepant and unparseable values are missing. Not available for GSE96058. Used only in the stage sensitivity feature sets, which are fit and evaluated on patients with recorded stage. |

---

## 3. Survival outcomes (harmonized)

The modeling pipeline reads survival from the harmonized `survival` table in
`data/processed/unified_cohorts.duckdb` (joined on `sample_id`). These are the
authoritative outcome variables used for all C-index/calibration/DCA analyses.

| Variable | Role | Units / encoding |
|---|---|---|
| `os_days` | OUTCOME (time) | Overall-survival follow-up time in **days** (continuous; DOUBLE). TCGA-BRCA uses TCGA-CDR `OS.time`. No minimum follow-up is applied; censored patients enter at their last follow-up. Rows with `os_days <= 0` are excluded from C-index calculations. |
| `os_event` | OUTCOME (event) | Overall-survival event indicator: **0 = censored, 1 = death** (BIGINT). |
| `dfs_days` | OUTCOME (time) | Disease-/recurrence-free survival time in **days** (DOUBLE). METABRIC RFS, GSE20685 relapse time where recorded, TCGA-BRCA TCGA-CDR progression-free interval; NULL for GSE96058. Used by the METABRIC recurrence sensitivity analysis. |
| `dfs_event` | OUTCOME (event) | Disease-/recurrence-free event indicator: **0 = censored, 1 = recurrence/relapse** (BIGINT). |
| `last_followup_days` | AUX | Last-contact time in days (DOUBLE), where available. |

Note on the legacy clinical CSVs: `data/legacy_v1_preprint/01_tcga_clinical.csv` and
`01_gse96058_clinical.csv` are the v1 (preprint) cohorts and are no longer read
by the pipeline. Their `time_to_event` column is in months.

---

## 4. PAM50-ROR comparator (baseline outputs)

Materialized in `data/processed/baselines_pam50.parquet` (long format, one row
per `cohort` x `sample_id` x `baseline`). Two baselines are present:
`PAM50_RORS_genefu` (genefu, `method_label = genefu_rorS`) and
`PAM50_ROR_surrogate` (`published_subtype_plus_proliferation`). The
**genefu ROR-S** baseline is the head-to-head comparator in the manuscript.

| Column | Role | Definition |
|---|---|---|
| `cohort` | ID | Cohort label (TCGA-BRCA, GSE96058, METABRIC, GSE20685). |
| `sample_id` | ID | Harmonized sample identifier (join key). |
| `baseline` | ID | `PAM50_RORS_genefu` or `PAM50_ROR_surrogate`. |
| `score` | BASELINE | Raw ROR/risk score on the method's native scale. |
| `score_0_100` | BASELINE | ROR score rescaled to **0-100** (verified min 0.0, max 100.0). This is the discrimination comparator score. |
| `method_label` | AUX | `genefu_rorS` (genefu research ROR-S) or `published_subtype_plus_proliferation` (surrogate). |
| `status` | AUX | `ok` (score computed) or `unavailable` (insufficient inputs). |
| `note` | AUX | Free-text provenance / unavailability reason. |
| `pam50_subtype_genefu` | BASELINE/AUX | genefu-called intrinsic subtype: one of {Basal, Her2, LumA, LumB, Normal}. |
| `rorS_risk` | BASELINE/AUX | genefu ROR-S risk group: {Low, Intermediate, High}. |
| `n_pam50_genes_input` | AUX | Number of PAM50 genes available to the scorer for that sample. |
| `er_status` | AUX | ER status used for risk-group cutpoints. |
| `tnbc_flag` | AUX | Triple-negative indicator (used to define the registered TNBC subgroup). |
| `intrinsic_subtype_pam50_published` | AUX | Published/source-provided intrinsic subtype, for cross-check. |

---

## 5. Per-file column tables

### 5a. `data/clinical/tcga_brca_clinical.csv` (TCGA-BRCA, written by `scripts/ingest_raw.py`)

| Column | Role | Notes |
|---|---|---|
| `sample_id`, `patient_id` | ID | TCGA primary-tumour barcode (-01) and patient barcode. |
| `age_at_dx` | AUX | `age_at_initial_pathologic_diagnosis`. |
| `stage_tnm` | AUX | `pathologic_stage`; mapped to `stage_ordinal`. |
| `er_status`, `pr_status`, `her2_status` | AUX | IHC receptor status; HER2 combines IHC and ISH. |
| `pam50_published` | AUX | `PAM50Call_RNAseq`. |
| `os_days`, `os_event` | OUTCOME | TCGA-CDR `OS.time` (days) and `OS`. |
| `pfi_days`, `pfi_event` | OUTCOME | TCGA-CDR progression-free interval. |

### 5b. `data/clinical/gse96058_clinical.csv` (GSE96058 / SCAN-B, written by `scripts/ingest_raw.py`)

| Column | Role | Notes |
|---|---|---|
| `sample_id`, `geo_accession` | ID | GEO sample title (F-number) and GSM accession. |
| `age_at_diagnosis` | AUX | Harmonized to `age_at_dx`. |
| `overall_survival_days`, `overall_survival_event` | OUTCOME | Parsed as `os_days`, `os_event`. |
| `er_status`, `pgr_status`, `her2_status`, `ki67_status`, `nhg` | AUX | Clinicopathology / receptor status. |
| `pam50_subtype` | AUX | Published PAM50 subtype. |
| `*_prediction_mgc`, `*_prediction_sgc` | AUX | SCAN-B multi-/single-gene classifier predictions. |
| `tumor_size`, `lymph_node_group`, `lymph_node_status` | AUX | Tumor burden / nodal status. |
| `endocrine_treated`, `chemo_treated` | AUX | Treatment flags. |

`data/clinical/legacy_v1_cohort_ids.csv` lists the sample IDs of the legacy
TCGA (213) and GSE96058 (1,483) sets; `qc.in_legacy_v1_cohort` flags them.

### 5c. `data/processed/04_pathway_features.parquet` (all cohorts; one row per scored sample)

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

See Section 4 for the full column table.

---

## 6. Harmonized `samples` and `survival` tables (`unified_cohorts.duckdb`)

The modeling frame in `scripts/train_ml_zoo.py` is built by joining
`04_pathway_features.parquet` to these two tables on `sample_id` / `cohort`.

**`samples`** (per sample): `sample_id` (ID), `cohort` (ID), `patient_id` (ID),
`age_at_dx` (SENSITIVITY INPUT), `sex`, `race`, `ethnicity`, `stage_tnm` (source of
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
| TCGA-BRCA | 1,064 | Training cohort (internal 5-fold CV; out-of-fold predictions used in pooled analyses) |
| GSE96058 / SCAN-B | 3,273 | External validation |
| METABRIC | 1,980 | External validation (also DFS/recurrence endpoint) |
| GSE20685 | 327 | External validation (OS only; no receptor status, cannot enter TNBC endpoint) |

Sample counts entering the genefu PAM50 ROR-S head-to-head and patient-
characteristics analyses differ from these scoring counts because of metric-
evaluability filters (non-positive survival time, model-specific missingness);
see `results/Table_1_denominator_crosswalk.csv`.
