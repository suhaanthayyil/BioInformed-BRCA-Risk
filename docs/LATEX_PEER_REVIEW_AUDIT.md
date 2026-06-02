# LaTeX Peer-Review Audit

Generated: 2026-06-02

Scope:

- Current repository head: `db054582ba5ee8995c88ce8c77aae60df6b1ffa7`
- Current branch state checked with `git status --short --branch`: clean and synced with `origin/main`.
- Current tracked LaTeX files: `results/Table_S1_clinical_baselines.tex` only.
- Historical full LaTeX manuscript source audited from commit `dbc8eda4c711784e97c3fc1fb715de3db6bfd0b8:paper/main.tex`.
- Historical bibliography audited from commit `dbc8eda4c711784e97c3fc1fb715de3db6bfd0b8:paper/main.bib`.
- Submitted Word manuscript previously audited separately from `/Users/suhaan/Downloads/manuscript (3).docx`.

Bottom line:

- The current GitHub repository no longer contains a full LaTeX manuscript source. The `paper/` and `manuscript/` directories were removed in commit `7cf98b5cb409ffb40fc59ff42a16710650ab2dd1`.
- The active tracked LaTeX file is only a supplementary clinical-baseline table fragment.
- The historical full LaTeX source is not submission-ready and does not match the submitted Word manuscript. It is a BMC template skeleton with TODO blocks.
- The central submitted result values in the repository pass `scripts/check_submission_consistency.py`.

## Current Repository LaTeX Inventory

Command evidence:

```text
find . -maxdepth 4 -type f \( -name '*.tex' -o -name '*.bib' -o -name '*.docx' \)
```

Result:

```text
./results/Table_S1_clinical_baselines.tex
```

Reviewer implication:

- If the paper or response letter says the repository contains manuscript LaTeX source, that is false for the current public tree.
- The repository should be described as code, results, figures, and traceability outputs, not as the active manuscript-source repository.

## Historical Full LaTeX Source Audit

Source audited:

```text
git show dbc8eda4c711784e97c3fc1fb715de3db6bfd0b8:paper/main.tex
```

### Lines 1-6

Status: structurally plausible but not current.

- Uses `\documentclass{bmcart}`, which is appropriate for a BMC-style LaTeX manuscript.
- Uses standard packages only.
- This file is not present in the current repository.

### Lines 13-19

Status: not submission-ready.

Issues:

- Line 13 title is different from the submitted Word manuscript title.
- Lines 15-19 contain TODO placeholders for email, affiliation, city, and country.
- This historical source cannot be treated as final metadata.

Reviewer risk:

- High if a reviewer compares this source to the submitted manuscript.

### Lines 22-36

Status: incomplete.

Issues:

- Abstract contains TODO blocks for Background, Results, and Conclusions.
- The Methods abstract is generic and does not include all core reported numbers.
- No citations are present, which is acceptable in an abstract for BMC guidance, but the rest of the manuscript also lacks citations.

### Lines 40-41

Status: incomplete.

Issue:

- Background section is only `[TODO: Suhaan writes Background.]`.

### Lines 45-46

Status: partially accurate, needs precision.

Text says:

```text
The harmonized database contained 4532 samples. The final Phase Three analysis set included 4003 samples with usable overall survival and official PAM50-ROR scores.
```

Repository-supported values:

- Harmonized total: `4,532`
- Patient-characteristics evaluable: `4,003`
- Official PAM50 head-to-head metric-evaluable: `4,002`
- Metric-evaluable events: `1,680`

Source:

- `results/Table_1_denominator_crosswalk.csv`
- `results/reports/harmonization_summary.json`

Recommended reviewer-proof wording:

```text
The harmonized database contained 4,532 samples. The patient-characteristics analysis set included 4,003 samples, while the official PAM50-ROR discrimination analysis included 4,002 metric-evaluable samples after excluding one METABRIC sample with non-positive survival time.
```

### Lines 48-49

Status: incomplete but directionally correct.

Issue:

- It says Figure 1 gives counts, but current repository no longer tracks `manuscript/figures/Figure_1.*`.
- In current GitHub, generated figure files are under `figures/`, not `manuscript/figures/`.

### Lines 54-55

Status: incomplete and under-specific.

Issues:

- It says seven pathway features were used, which is correct.
- It does not specify the exact pathway component count.
- Current pathway-scoring metadata says:
  - Hallmark gene sets: `50`
  - Reactome curated: `244`
  - KEGG MEDICUS curated: `50`
  - Mandatory aggregation sets: `29`
  - Total gene sets scored: `347`

Source:

- `data/processed/pathway_scoring.metadata.json`

Recommended wording:

```text
Pathway scoring used MSigDB v2024.1.Hs. Across cohorts, 347 gene sets were scored, including 50 Hallmark, 244 curated Reactome, and 50 KEGG MEDICUS sets; 29 pre-specified component sets were aggregated into seven pathway features.
```

### Lines 57-58

Status: contains a substantive mismatch.

Issue:

- Historical LaTeX says surrogate Oncotype DX and MammaPrint scores were retained as supplementary comparators.
- Current repository indicates MammaPrint, EndoPredict, and Breast Cancer Index were not computed because exact public formula/coefficient sets were unavailable.

Source:

- `results/reports/clinical_baselines_qc.md`
- `results/Table_S1_clinical_baselines.csv`

Recommended wording:

```text
Official PAM50 subtype and ROR-S risk were computed using genefu in R. A surrogate Oncotype DX 21-gene score was retained as a supplementary comparator where ER status and gene coverage allowed. MammaPrint, EndoPredict, and Breast Cancer Index were not computed because validated public coefficient sets were not available from the public inputs, and surrogate scores were not fabricated.
```

### Lines 60-61

Status: high-level accurate but too sparse for reviewer reproduction.

Issue:

- The model list is correct at a high level.
- It does not report the actual selected model configurations.

Repository-supported selected configs:

- Cox PH: `penalizer = 0.01`
- Elastic Net Cox: `alpha = 0.001`, `l1_ratio = 0.1`
- Random Survival Forest: `n_estimators = 500`, `max_depth = 4`, `min_samples_leaf = 1`, `max_features = sqrt`
- Gradient Boosted Survival: `n_estimators = 100`, `learning_rate = 0.03`, `max_depth = 1`, `subsample = 1.0`
- DeepSurv: hidden dimensions `[128, 64]`, dropout `0.2`, learning rate `0.0005`, weight decay `0.0001`, epochs `150`, patience `20`
- Stacked ensemble: `Gradient_Boosted_Survival`, `DeepSurv`, `Elastic_Net_Cox`

Source:

- `results/Table_S2_ml_internal_cv.csv`
- `scripts/train_ml_zoo.py`

Reviewer risk:

- High if the manuscript gives different hyperparameters from these values.

### Lines 63-64

Status: numerically correct.

Text:

```text
The headline model had internal C-index 0.642 and delta versus Cox +0.0420.
```

Repository evidence:

- Exact headline internal CV C-index: `0.6418257777973452`
- Exact delta versus Cox: `0.041963998078774156`

Source:

- `results/Table_S2_ml_internal_cv.csv`

Verdict:

- Correct after rounding.

### Lines 66-67

Status: correct at high level.

Notes:

- External cohorts are GSE96058/SCAN-B, METABRIC, and GSE20685.
- TCGA-BRCA was the training cohort.
- This matches repository implementation.

Source:

- `results/Table_2_external_validation.csv`
- `scripts/external_validation.py`

### Lines 69-70

Status: mostly correct but under-cited.

Issues:

- Statistical method list is consistent with repository output tables.
- The phrase `decile-smoothed integrated calibration index` is not clearly tied to the exact implementation and should be checked against `scripts/calibration_analysis.py`.
- Bibliographic support exists for DerSimonian-Laird, ICI, DCA, and NRI/IDI, but historical LaTeX has no inline `\cite{}` commands.

### Lines 72-73

Status: central numbers correct, commit SHA stale relative to current head.

Correct values:

- Primary endpoint: TNBC meta-analyzed delta C-index
- Success threshold: `delta >= +0.03`, `p < 0.05`
- Observed delta: `+0.0144`
- 95% CI: `[-0.0472, +0.0760]`
- p-value: `0.6466`
- Endpoint status: not met

Commit note:

- Line 73 uses commit `d8b5ee75ee9f091b939208f8ddc811ca7e048527`, which exists and corresponds to the same timestamp/message as `7783f4f0bfaaa6bdc611c78d33ccda621c6b243d`.
- Current README uses `7783f4f0bfaaa6bdc611c78d33ccda621c6b243d`.
- Use one canonical SHA consistently, preferably `7783f4f0bfaaa6bdc611c78d33ccda621c6b243d`.

### Lines 75-76

Status: stale.

Issues:

- Repository URL placeholder remains.
- Commit SHA `a704e9ff339ef765f9430baf859a27d0a4b2e001` is not current.
- Current head is `db054582ba5ee8995c88ce8c77aae60df6b1ffa7`.

Recommended wording:

```text
All code and derived outputs are available at https://github.com/suhaanthayyil/BioInformed-BRCA-Risk. The repository state aligned to the submitted manuscript is commit db054582ba5ee8995c88ce8c77aae60df6b1ffa7.
```

### Lines 81-97

Status: incomplete.

Issues:

- Results, Discussion, and Conclusions are TODO blocks.
- This historical LaTeX source is not a complete manuscript.

### Lines 99-120

Status: partially complete but TODO-laden.

Issues:

- Declarations are present, which aligns with BMC requirements.
- Competing interests, funding, author contributions, acknowledgements, and repository URL contain TODO placeholders.
- Current submitted Word manuscript has fuller declarations than this historical LaTeX.

### Lines 122-123

Status: bibliography is not wired into text.

Issue:

- `\bibliography{main}` exists.
- The historical `.tex` contains no `\cite{...}` commands.
- Compiling this source would not produce meaningful in-text numbered citations.

Reviewer risk:

- High if this LaTeX is submitted or shared as manuscript source.

## Historical Bibliography Audit

Source audited:

```text
git show dbc8eda4c711784e97c3fc1fb715de3db6bfd0b8:paper/main.bib
```

Structural findings:

- `18` BibTeX entries.
- No DOI fields.
- No PMID fields.
- No URL fields.
- No `TODO_BIB` markers.
- All entries are unused by `paper/main.tex` because the `.tex` has no `\cite{}` commands.

Citation metadata issues:

- `pereira2016somatic`: title uses `refines`; the corrected Nature Communications title uses `refine`.
- `tcga2012comprehensive`: DOI should be added: `10.1038/nature11412`.
- `paik2004gene`: DOI should be added: `10.1056/NEJMoa041588`.
- `cox1972regression`: DOI should be added: `10.1111/j.2517-6161.1972.tb00899.x`.
- `ishwaran2008random`: DOI should be added: `10.1214/08-AOAS169`.
- `derSimonian1986meta`: DOI should be added: `10.1016/0197-2456(86)90046-2`.
- `hothorn2006survival`: DOI should be added: `10.1093/biostatistics/kxj011`.
- `austin2019graphical`: DOI should be added: `10.1002/sim.8570`.
- `collins2015tripod`: if citing Annals of Internal Medicine TRIPOD statement, DOI should be `10.7326/M14-0697`. Crossref may also return the parallel Journal of Clinical Epidemiology version; choose the one matching the journal listed.

Verdict:

- Bibliography is salvageable but not reviewer-proof without DOI fields and inline citations.

## Current Active LaTeX Table Audit

Source:

```text
results/Table_S1_clinical_baselines.tex
```

### Lines 1-4

Status: fragment only.

Issues:

- This is a raw `tabular`, not a standalone table environment.
- No caption, label, resize, or `sidewaystable`.
- Column count is very wide: `llllrrrrrrrrrrrrrrll`.
- This will not fit normal manuscript page width without resizing or rotating.

### Lines 3-40

Status: numeric values match the CSV, but LaTeX escaping is broken.

Issue:

- Raw underscores appear throughout identifiers:
  - `Oncotype_DX_21_gene_surrogate`
  - `PAM50_ROR_official`
  - `PAM50_ROR_surrogate`
  - `non_tnbc`
  - `too_few_events`

LaTeX risk:

- These will cause compile errors in normal text mode unless escaped as `\_` or wrapped in `\texttt{}`.

### Lines 5-13

Status: numerically correct.

Notes:

- GSE20685 Oncotype and PAM50 surrogate are not evaluable.
- GSE20685 official PAM50-ROR overall C-index `0.635` matches the CSV rounded value.

### Lines 14-22

Status: numerically correct, wording partly stale.

Notes:

- GSE96058 Oncotype and PAM50 rows match `results/Table_S1_clinical_baselines.csv`.
- Surrogate PAM50 note says official genefu should replace the surrogate if package scoring succeeds. Official genefu scoring has already succeeded, so this note should be rewritten.

### Lines 23-31

Status: numerically correct, wording partly stale.

Notes:

- METABRIC official PAM50-ROR row is correct in this table: n `1980`, events `1143`, Harrell C `0.595`.
- Surrogate PAM50 note is stale for the same reason above.

### Lines 32-40

Status: numerically correct, wording partly stale.

Notes:

- TCGA-BRCA official PAM50-ROR row is correct.
- Surrogate PAM50 note is stale.

## Current Repo Alignment Audit

Current result verification:

```text
python3 scripts/check_submission_consistency.py
```

All central checks pass:

- TNBC primary delta: `0.014398`
- TNBC primary p-value: `0.646627`
- Primary endpoint status: not met
- Headline internal CV C-index: `0.641826`
- Headline internal delta versus Cox: `0.041964`
- External Harrell C values for GBSA and PAM50-ROR across TCGA-BRCA, GSE96058, METABRIC, and GSE20685
- Mean 5-year ICI values
- Mean DCA net-benefit delta
- Stability mean Spearman rho

Current traceability files:

- `docs/SUBMITTED_MANUSCRIPT_REPOSITORY_CROSSWALK.md`
- `results/RESULTS_LOG.md`
- `results/Table_1_denominator_crosswalk.csv`
- `results/Table_S2_internal_cv_manuscript_crosswalk.csv`
- `scripts/check_submission_consistency.py`

Reviewer-proofing issue:

- `README.md` says `Analysis set: 4,003 patients with usable overall survival and official PAM50-ROR scores`. For precision, update this to distinguish:
  - `4,003` patient-characteristics evaluable
  - `4,002` metric-evaluable for official PAM50-ROR discrimination

## Wording Review

Safe wording:

- `competitive with PAM50-ROR`
- `comparable to PAM50-ROR`
- `did not significantly exceed PAM50-ROR`
- `primary endpoint was not met`
- `exploratory secondary analysis`
- `transparent sensitivity analysis`

Risky wording:

- `superior to PAM50-ROR`, unless describing the failed pre-registered hypothesis.
- `outperforms PAM50`, unless explicitly negated.
- `pre-specified 391-variant sensitivity analysis`, unless a separate timestamped pre-specification file exists.
- `all consistent in direction across cohorts` for calibration, DCA, and IDI together. IDI is positive in each cohort, but DCA and ICI have cohort-level exceptions.

## Reviewer-Response Ready Counts

Use these exact denominator statements:

- Harmonized database: `4,532` samples.
- Patient-characteristics analysis set: `4,003` samples.
- Official PAM50-ROR discrimination metric set: `4,002` samples.
- Metric-evaluable event count: `1,680`.
- TCGA-BRCA: `213` samples, `132` events.
- GSE96058: `1,483` samples, `322` events.
- METABRIC: `2,509` harmonized, `1,980` patient-characteristics evaluable, `1,979` metric-evaluable, `1,144` patient-characteristics events, `1,143` metric-evaluable events.
- GSE20685: `327` samples, `83` events.
- TNBC primary endpoint: `408` patients, `214` events.
- External-only TNBC sensitivity: `374` patients, `195` events.

## Final Audit Verdict

The repository's central numerical outputs are internally consistent and support the submitted negative primary endpoint framing.

However, the LaTeX situation is not reviewer-proof:

1. There is no current full manuscript LaTeX source in GitHub.
2. The historical full LaTeX source is incomplete and stale.
3. The historical bibliography has no DOI fields and no inline citations.
4. The only current tracked LaTeX table fragment has compile-breaking raw underscores and stale surrogate-PAM50 notes.

Recommended action:

- Treat the current GitHub repo as a code/results repository, not a manuscript-source repository.
- If LaTeX source is needed, recreate it from the submitted Word manuscript rather than using the historical `paper/main.tex`.
- Fix or remove `results/Table_S1_clinical_baselines.tex` before pointing reviewers to TeX tables.
