# Registered plan versus reported analysis

Registration: `docs/PRIMARY_ENDPOINT.md`, commit `7783f4f` (2026-05-16). Machine-readable copy: `results/Table_S20_registration_deviations.csv`.

## Cohorts

- Registered: At least three independent cohorts (not named)
- Submitted manuscript: Submitted: TCGA-BRCA (training), GSE96058, METABRIC, GSE20685
- Revised analysis: Revised: same four cohorts
- Reason: METABRIC and GSE20685 were obtained after registration. TCGA-BRCA and GSE96058 had already been analysed in the April 2026 preprint, so only METABRIC and GSE20685 are validation cohorts that were unseen at registration.. The registration did not name the training cohort or the PAM50 implementation; TCGA-BRCA as sole training cohort and genefu::rorS were fixed after registration

## Cohort inclusion

- Registered: Not specified
- Submitted manuscript: Submitted: v1 cohorts reused (TCGA n = 213; GSE96058 n = 1,483). These excluded patients censored before 5 years, and included 37 TCGA normal-tissue and 3 metastatic samples and 60 GSE96058 technical replicates
- Revised analysis: Revised: all female primary-tumour TCGA-BRCA patients with TCGA-CDR OS (n = 1,064) and all non-replicate GSE96058 patients with OS (n = 3,273); no follow-up restriction
- Reason: The v1 cohorts were built for a binary 5-year outcome and were carried into the survival analysis by mistake; corrected after the journal's pre-publication check

## Model inputs

- Registered: Seven pathway scores
- Submitted manuscript: Submitted: seven pathways + age + ordinal stage; stage was missing (median-imputed) for every external patient because of a coding error, and TCGA 'No_Conversion' was coded as stage I
- Revised analysis: Revised: seven pathway scores only, as registered
- Reason: Restores the registered feature set. The endpoint is all-cause overall survival and PAM50 ROR-S contains no age term, so adding age to the ML model would confound the comparison. Pathways + age and pathways + age + correctly coded stage (fit and evaluated only where stage is recorded; GSE96058 releases none) are reported as sensitivity analyses (Table 6, Table S5)

## Pathway definitions

- Registered: Seven interpretable biological pathway scores
- Submitted manuscript: Seven MSigDB v2024.1 aggregates: immune, proliferation, DNA repair, metabolism, stromal/EMT, apoptosis/stress, hormone
- Revised analysis: Unchanged from submission
- Reason: Defined after registration to use versioned public gene sets; the preprint's 45-gene sets (including HER2 and angiogenesis axes) were replaced

## Candidate models

- Registered: RSF, GBSA or DeepSurv, chosen by internal CV
- Submitted manuscript: Six models including Cox PH, elastic-net Cox and a stacked ensemble; headline = best internal CV over all six (GBSA)
- Revised analysis: Same six models; headline restricted to the registered candidates (RSF, GBSA, DeepSurv) and chosen by internal CV; linear models and the ensemble are comparators
- Reason: Restores the registered selection rule; linear baselines are needed for secondary claim (a)

## Internal cross-validation

- Registered: Internal CV on the training cohort
- Submitted manuscript: DeepSurv early stopping monitored the held-out test fold, inflating its CV C-index; DeepSurv ran on the Apple MPS device
- Revised analysis: DeepSurv early stopping uses a 20% split of the training fold only; CPU, single-threaded, fixed seed
- Reason: Removes test-fold leakage and makes DeepSurv reproducible

## Hyperparameter search

- Registered: Winner by internal CV (grid not specified)
- Submitted manuscript: Compact local grid (RSF 18, GBSA 16, DeepSurv 4 configurations), recorded as not the full protocol grid
- Revised analysis: Expanded grid (RSF 75, GBSA 81, DeepSurv 40 configurations); top 5 per model re-scored with 5x5 repeated CV; selection on TCGA only; external cohorts scored once
- Reason: The compact grid was flagged in the metadata as a local shortcut; a wider search and repeated CV give each registered model a fair tuning budget and reduce split-to-split noise (150 training deaths). Four revision rounds (0-3) were run, including a pathways + age round that was not used (docs/HYPERPARAMETERS.md); selection always used TCGA internal CV only, but the external cohorts were not blind during the revision because their earlier results had been seen. External results were similar across rounds; round 3 is final

## PAM50 comparator

- Registered: PAM50 / PAM50-ROR
- Submitted manuscript: genefu::rorS described as 'official PAM50-ROR'
- Revised analysis: genefu::rorS, a research implementation of PAM50 ROR-S (subtype correlations + proliferation), not the Prosigna assay
- Reason: The clinical Prosigna test cannot be computed from public expression data; terminology corrected

## Oncotype DX comparison

- Registered: Non-inferior overall and superior in TNBC versus Oncotype DX surrogate
- Submitted manuscript: Surrogate scored but only reported as a standalone baseline
- Revised analysis: Head-to-head versus the 21-gene research surrogate in ER-positive patients (Table S18); TNBC comparison reported as not evaluable
- Reason: The 21-gene recurrence score is defined for ER-positive disease, so the comparison population changes from the overall cohort (as registered) to ER-positive patients, and TNBC superiority cannot be tested

## Primary endpoint: training cohort

- Registered: TNBC delta C-index meta-analysed across all cohorts with TNBC patients
- Submitted manuscript: TCGA-BRCA contributed in-sample predictions from the final fitted model
- Revised analysis: TCGA-BRCA contributes out-of-fold cross-validated predictions; an external-only meta-analysis is reported alongside
- Reason: In-sample predictions on the training cohort are optimistic

## Secondary (a) ML vs Cox >= +0.02

- Registered: Registered
- Submitted manuscript: Internal CV only
- Revised analysis: Internal CV and each external cohort (Table S18)
- Reason: Reported in full

## Secondary (b) immune pathway dominant in TNBC and HER2+

- Registered: Registered
- Submitted manuscript: Not reported
- Revised analysis: Permutation importance by subgroup in external cohorts (Table S18)
- Reason: Previously omitted

## Secondary (c) combined ML + PAM50 + clinical IDI > 0.05

- Registered: Registered
- Submitted manuscript: IDI of ML vs PAM50 reported instead
- Revised analysis: Combined (ML + PAM50 ROR-S + age) vs PAM50 ROR-S IDI at 5 years (Table S18)
- Reason: Previously the comparison was swapped

## Secondary (d) non-inferiority to PAM50, margin -0.02

- Registered: Registered
- Submitted manuscript: Not tested
- Revised analysis: One-sided 95% lower bound of the pooled external delta (Table S18)
- Reason: Previously omitted

## Post-hoc analyses

- Registered: None registered
- Submitted manuscript: Rescue model search, calibration, DCA, NRI/IDI, within-subtype, stability, ablation, scoring sensitivity, recurrence endpoint, learning curve
- Revised analysis: Same, plus locked-model sensitivity on full cohorts (Table S19)
- Reason: All exploratory and labelled as such. The rescue search (391 candidates selected on TCGA CV) chose a Cox model on the seven pathways plus age and stage; its TNBC advantage over PAM50 ROR-S (external-only delta +0.13) reflects the clinical covariates, consistent with Table S5, and is not evidence for the pathway model

## Registration document

- Registered: 'Modifications after analysis begins: NOT PERMITTED'
- Submitted manuscript: File edited 2026-05-24 and 2026-06-22; its timeline note wrongly said no results preceded it
- Revised analysis: Timeline corrected; this table lists every deviation
- Reason: Transparency
