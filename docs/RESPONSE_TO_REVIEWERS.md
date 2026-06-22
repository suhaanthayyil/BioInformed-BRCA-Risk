# Response to Reviewers — BMC Artificial Intelligence (Major Revision)

**Manuscript:** *Pathway-based machine learning is competitive with but does
not exceed PAM50-ROR for breast cancer prognosis: a pre-registered multi-cohort
benchmark in 4,532 patients.*

This is a point-by-point scaffold. For each editor (E) and reviewer (R)
comment: a one-line restatement of the ask, a **Repo action** line pointing at
the repository artifact that addresses it, and a `[PROSE TODO — author]`
placeholder where the manuscript text is the author's job.

The headline finding is unchanged and reported transparently: the pre-registered
TNBC superiority endpoint was **not met** (meta delta C-index +0.0144, 95% CI
[-0.0472, +0.0760], p = 0.6466); the conclusion is "competitive, not superior."
See `docs/PRIMARY_ENDPOINT.md`, `docs/STORY.md`, and
`data/processed/phase6_primary_endpoint.json`.

---

## Editor comments

### E1 — Clinician / biostatistician validation of outcomes and secondary metrics
**Ask:** Discuss involving clinicians and biostatisticians to validate the
outcome definitions and secondary metrics.
**Repo action:** Paper-only; no repository artifact. The endpoint definitions
and metric code are traceable via `docs/DATA_DICTIONARY.md` and
`docs/PRIMARY_ENDPOINT.md` for any such review.
`[PROSE TODO — author]`

### E2 — Broader recent-literature review
**Ask:** Expand the intro/discussion with a broader, more recent literature
review.
**Repo action:** Add citations to `paper/main.bib` (including the two R2.4 DOIs,
see R2.4 below) and expand intro/discussion.
`[PROSE TODO — author]`

### E3 — Source code + trained weights + curated testing datasets (paramount)
**Ask:** Provide source code, trained model weights, and curated testing
datasets so the work is reproducible.
**Repo action:**
- Trained weights now committed under `models/` with a `models/README.md`
  describing each artifact (manuscript model name, training cohort TCGA-BRCA,
  feature list, seed 42, CV C-index, library versions, per-file SHA256).
- Pipeline intermediates committed under `data/processed/`.
- `scripts/predict.py` loads the committed weights cold (importable wrappers)
  and emits risk scores from a feature matrix.
- `tests/test_smoke.py` recomputes a headline number (not a CSV diff).
- Reproducibility infra added: `Makefile`, `Dockerfile`, and
  `.github/workflows/ci.yml`.
`[PROSE TODO — author]` (Code/Data Availability statement with clone URL + DOI.)

---

## Reviewer 1

### R1.1 — Frame calibration/DCA/IDI/NRI as secondary/exploratory
**Ask:** Present calibration, DCA, IDI, and NRI as secondary/exploratory; do not
imply demonstrated clinical utility.
**Repo action:** `analysis_type='exploratory_secondary'` columns added to
`results/Table_S10_calibration.csv`, `results/Table_S11_dca.csv`, and
`results/Table_S8_within_subtype_external.csv`; figure captions to be marked
exploratory.
`[PROSE TODO — author]` (mark `fig:dca` and `fig:calibration` captions
exploratory.)

### R1.2 — Small training cohort: stability, feature-rank, generalizability; justify TCGA-only training
**Ask:** Address the small TCGA-only training cohort — stability, feature-rank
robustness, generalizability — and justify training on TCGA only.
**Repo action:** `scripts/learning_curve.py` +
`results/Table_S17_learning_curve.csv` + `figures/fig_learning_curve.*`
(training-size sensitivity, plus a train-on-METABRIC-validate-on-rest
cross-check). Existing stability evidence in `results/Table_S13_stability.csv`.
`[PROSE TODO — author]`

### R1.3 — Distinguish "no demonstrated superiority" from "no effect"; acknowledge CI width
**Ask:** Do not equate failure to demonstrate superiority with absence of
effect; acknowledge the width of the confidence interval.
**Repo action:** Forest plots now annotated with event counts; the primary CI
([-0.0472, +0.0760]) is reported alongside the point estimate. Lack of
superiority is not no-effect.
`[PROSE TODO — author]`

### R1.4 — Expand harmonization detail
**Ask:** Expand harmonization: probe mapping, cross-cohort gene overlap,
platform handling, residual heterogeneity.
**Repo action:** `docs/HARMONIZATION.md` +
`results/reports/cohort_gene_overlap.csv`; expanded per-cohort download scripts.
Variable definitions in `docs/DATA_DICTIONARY.md`; source terms in
`docs/DATA_USE.md`.
`[PROSE TODO — author]`

### R1.5 — Would the metric improvements change a clinical decision?
**Ask:** Would the observed metric improvements actually change a clinical
decision?
**Repo action:** DCA/calibration analyses surface the decision threshold
(0.20 5-year mortality) and absolute net benefit (mean 5-year net-benefit delta
vs PAM50-ROR = +0.0156 across cohorts/thresholds; `results/Table_S11_dca.csv`).
`[PROSE TODO — author]` (clinical interpretation of absolute net benefit.)

### R1.6 — Hormone / age / DNA-repair are descriptive, not mechanistic
**Ask:** Treat the top features (Hormone, age, DNA-repair) as descriptive
signals, not mechanistic findings.
**Repo action:** Stale feature matrix fixed so `Pathway_Hormone` and
`Pathway_DNA_Repair` are reproducible from the committed
`data/processed/04_pathway_features.parquet`; stability in
`results/Table_S13_stability.csv` (top-3 by frequency: Pathway_Hormone 92%,
age_at_dx 75%, Pathway_DNA_Repair 62%).
`[PROSE TODO — author]` (descriptive-not-mechanistic framing.)

### R1-minor — Figure legends, retrospective/public-data caveat, terminology, prognosis ≠ prediction
**Ask:** Figure legends should state train/val/pooled, n, and events; add a
retrospective/public-data caveat; standardize terminology (PAM50 vs PAM50-ROR
vs genefu); clarify prognosis is not treatment-response prediction.
**Repo action:** Figures annotated with event counts and train/val/pooled
labels; data-use/de-identification note in `docs/DATA_USE.md`.
`[PROSE TODO — author]` (terminology pass; prognosis-not-prediction sentence.)

---

## Reviewer 2

### R2.1 — PAM50-ROR is a recurrence tool; benchmark is OS — stress-test this
**Ask:** PAM50-ROR was developed for recurrence; the benchmark uses overall
survival. Stress-test the endpoint mismatch.
**Repo action:** `scripts/recurrence_sensitivity.py` +
`results/Table_S16_metabric_recurrence_sensitivity.csv` (locked GBSA vs official
PAM50-ROR on the METABRIC disease-/recurrence-free endpoint, using the
harmonized `dfs_days`/`dfs_event` available for METABRIC).
`[PROSE TODO — author]`

### R2.2 — Generalizability of TCGA-only training
**Ask:** Demonstrate generalizability given TCGA-only training.
**Repo action:** `scripts/learning_curve.py` +
`results/Table_S17_learning_curve.csv` + `figures/fig_learning_curve.*`,
including the train-on-METABRIC cross-check.
`[PROSE TODO — author]`

### R2.3 — More hyperparameter-selection detail
**Ask:** Provide more detail on hyperparameter search and model selection.
**Repo action:** `docs/HYPERPARAMETERS.md` — per-model search spaces (from
`scripts/train_ml_zoo.py` `model_spaces()`), the 5-fold internal stratified CV
selection rule on TCGA-BRCA maximizing Harrell C-index (seed 42), the selected
configs (from `results/Table_S2_ml_internal_cv.csv`), and the recorded
`grid_mode = compact_local` caveat.
`[PROSE TODO — author]` (none strictly required; cite `docs/HYPERPARAMETERS.md`.)

### R2.4 — Discuss two specific references
**Ask:** Discuss `10.1186/s12885-024-12331-5` and `10.1038/s41598-025-21746-4`.
**Repo action:** Add both DOIs to `paper/main.bib` and cite them in the
intro/discussion. Candidate BibTeX entries below — **verify all fields
(authors, title, year, volume, pages) against the publisher record before
use.**

```bibtex
@article{r2_4_bmccancer_2024,
  title   = {VERIFY: title for DOI 10.1186/s12885-024-12331-5},
  author  = {VERIFY: author list},
  journal = {BMC Cancer},
  year    = {2024},
  doi     = {10.1186/s12885-024-12331-5},
  note    = {Verify full bibliographic fields against the publisher record.}
}

@article{r2_4_scirep_2025,
  title   = {VERIFY: title for DOI 10.1038/s41598-025-21746-4},
  author  = {VERIFY: author list},
  journal = {Scientific Reports},
  year    = {2025},
  doi     = {10.1038/s41598-025-21746-4},
  note    = {Verify full bibliographic fields against the publisher record.}
}
```
`[PROSE TODO — author]` (discuss both works; confirm BibTeX fields.)

### R2.5 — Typo "monotonically" → "consistently"
**Ask:** Fix the word "monotonically" (also factually loose for a 3-of-4 count).
**Repo action:** Edit at `paper/main.tex:200`.
`[PROSE TODO — author]` (apply the single-word edit and recompile.)

---

## Reviewer 3

### R3.1 — PAM50-ROR is a recurrence tool; OS benchmark mismatch (same as R2.1)
**Ask:** Same endpoint-mismatch concern as R2.1.
**Repo action:** `scripts/recurrence_sensitivity.py` +
`results/Table_S16_metabric_recurrence_sensitivity.csv`.
`[PROSE TODO — author]`

### R3.2 — Generalizability / small training cohort (same as R1.2 / R2.2)
**Ask:** Generalizability given the small TCGA training set.
**Repo action:** `scripts/learning_curve.py` +
`results/Table_S17_learning_curve.csv` + `figures/fig_learning_curve.*`.
`[PROSE TODO — author]`

### R3.3 — Exploratory framing of secondary metrics (same as R1.1)
**Ask:** Frame calibration/DCA/IDI/NRI as exploratory.
**Repo action:** `analysis_type='exploratory_secondary'` columns on
`results/Table_S10_calibration.csv`, `results/Table_S11_dca.csv`,
`results/Table_S8_within_subtype_external.csv`; captions marked exploratory.
`[PROSE TODO — author]`

### R3.4 — Descriptive, not mechanistic, feature interpretation (same as R1.6)
**Ask:** Treat top features as descriptive, not mechanistic.
**Repo action:** Stale feature matrix fixed so `Pathway_Hormone` /
`Pathway_DNA_Repair` are reproducible; stability in
`results/Table_S13_stability.csv`.
`[PROSE TODO — author]` (descriptive-not-mechanistic framing.)
