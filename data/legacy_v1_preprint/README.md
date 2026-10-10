# Legacy v1 cohorts (April 2026 preprint). Not used by the current analysis.

These files define the cohorts of the original two-cohort, binary-outcome
analysis (notebooks 01-14, preprint April 2026). They are kept only so the
preprint can be traced.

Selection rule used then: a patient was kept only if they died within 5 years
or were followed for at least 5 years. Patients alive and censored before
5 years were excluded, which is why no living patient here has less than
60 months of follow-up.

Further problems with these sets:
- `01_tcga_clinical.csv` / `02_tcga_feature_matrix.csv` (213 rows) contain 173
  primary tumours, 37 solid-tissue-normal samples (barcode -11) and 3 metastatic
  samples (-06); 40 patients appear twice.
- `01_gse96058_clinical.csv` (1,483 rows) contains 60 technical replicates
  (sample names ending in "repl").
- `ablation_results.csv` is the v1 classifier ablation.

The current cohorts are built by `scripts/ingest_raw.py` from the raw sources
with no follow-up restriction: `data/clinical/tcga_brca_clinical.csv`
(1,064 patients) and `data/clinical/gse96058_clinical.csv` (3,273 patients).
