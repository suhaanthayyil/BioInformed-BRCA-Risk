# Locked v1 artifacts (2026-05-16). Not used by the current analysis except Table S19.

These are the model weights and inputs exactly as locked on 2026-05-16 and used in
the originally submitted manuscript. Their SHA256 digests match the original
manifest (see `docs/MODEL_SHA256.txt`).

- Training cohort: the legacy TCGA-BRCA set of 213 samples (see
  `data/legacy_v1_preprint/README.md`). It kept only patients who died within
  5 years or were followed for at least 5 years, so patients censored before
  5 years were excluded; it also contained 37 solid-tissue-normal and 3
  metastatic samples.
- `unified_cohorts_v1.duckdb` is the harmonized database of that version and
  shows the same pattern for TCGA-BRCA and GSE96058 (no living patient followed
  for less than about 60 months).
- Features: seven pathways + age + stage; stage was missing (median-imputed) for
  every external patient because of a coding error.

`scripts/locked_model_sensitivity.py` scores these weights on every eligible
patient in the current cohorts (Table S19). Nothing else reads this folder.
