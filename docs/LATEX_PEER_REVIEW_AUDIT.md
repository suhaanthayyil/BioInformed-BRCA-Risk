# LaTeX Peer-Review Audit

Generated: 2026-06-02

Current repository state after remediation:

- Full manuscript LaTeX source: `paper/main.tex`
- DOI-backed BibTeX bibliography: `paper/main.bib`
- Compiled manuscript PDF: `paper/main.pdf`
- Word manuscript export: `paper/main.docx`
- Active clinical-baseline LaTeX table fragment: `results/Table_S1_clinical_baselines.tex`
- Submitted-value consistency checker: `scripts/check_submission_consistency.py`

## Resolution Summary

Earlier audit finding: "Current repo cannot honestly be described as containing manuscript LaTeX."

Current status: resolved. The repository now contains `paper/main.tex` and `paper/main.bib`.

Earlier audit finding: historical `paper/main.tex` was incomplete, stale, and full of TODOs.

Current status: resolved for the active repository. `paper/main.tex` was recreated from the supplied manuscript LaTeX and cleaned for reviewer-facing traceability. It has no `TODO` placeholders.

Earlier audit finding: historical bibliography had no inline `\cite{}` usage and no DOI fields.

Current status: resolved. `paper/main.tex` uses BibTeX `\cite{...}` calls, and `paper/main.bib` contains DOI fields plus LaTeX-printable DOI URL notes for all 28 cited references.

Earlier audit finding: `results/Table_S1_clinical_baselines.tex` contained raw underscores that could break LaTeX compilation.

Current status: resolved. The table fragment was regenerated with LaTeX escaping, so identifiers such as `PAM50_ROR_official` are emitted as `PAM50\_ROR\_official`.

Earlier audit finding: Table S1 contained stale notes saying official genefu should replace surrogate PAM50.

Current status: resolved. `results/Table_S1_clinical_baselines.csv` and `.tex` now state that official genefu PAM50-ROR is the primary comparator and the surrogate is retained only for historical comparison.

Earlier audit finding: core numbers pass the consistency script.

Current status: still true. `python3 scripts/check_submission_consistency.py` passes all central values.

## Current LaTeX Source Checks

Checked file: `paper/main.tex`

Static checks:

- No `TODO` placeholders.
- No remaining manual numeric citation pattern of the form `~[1]` or `~[1, 2]`.
- All `\cite{...}` keys are present in `paper/main.bib`.
- Figure paths point to existing files in `figures/` through `\graphicspath{{../figures/}{figures/}}`.
- The Methods model-configuration paragraph matches `results/Table_S2_ml_internal_cv.csv`.
- The primary endpoint remains not met and is not reframed as superiority.

Compile checks:

- `tectonic --keep-logs main.tex` completed successfully from `paper/`.
- `paper/main.pdf` was produced with resolved citations and no LaTeX errors.
- Remaining TeX warnings are layout warnings from long table/pathway labels, not missing citations, missing figures, or fatal compilation problems.
- The compiled PDF bibliography visibly prints DOI URLs.
- `python3 scripts/build_manuscript_docx.py` produced `paper/main.docx` from `paper/main.tex` and `paper/main.bib`.
- The DOCX was rendered through LibreOffice using `render_docx.py`; the title block, table-heavy pages, and reference pages were visually spot-checked.

## Corrected Model Configuration Text

The fixed LaTeX reports the repository-selected model configurations:

- Cox PH: L2 penalizer `0.01`
- Elastic-net Cox: alpha `0.001`, l1-ratio `0.1`
- Random Survival Forest: 500 trees, max depth `4`, min leaf `1`, max features `sqrt`
- Gradient-Boosted Survival Analysis: 100 estimators, learning rate `0.03`, max depth `1`, subsample `1.0`
- DeepSurv: hidden dimensions `[128, 64]`, dropout `0.2`, learning rate `0.0005`, weight decay `0.0001`, 150 epochs, patience 20
- Stacked ensemble: Gradient-Boosted Survival Analysis, DeepSurv, and elastic-net Cox with a Cox meta-learner

Source: `results/Table_S2_ml_internal_cv.csv` and `scripts/train_ml_zoo.py`.

## Corrected Denominator Language

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

Source: `results/Table_1_denominator_crosswalk.csv`, `results/Table_1_patient_characteristics.csv`, `results/Table_2_external_validation.csv`, and `results/Table_3_head_to_head.csv`.

## Primary Endpoint Integrity

The fixed LaTeX keeps the locked primary endpoint unchanged:

- Model: Gradient-Boosted Survival Analysis
- Comparator: official `genefu::rorS` PAM50-ROR
- Endpoint: TNBC random-effects delta C-index
- Threshold: delta C-index `>= +0.03`, `p < 0.05`
- Result: delta `+0.0144`
- 95% CI: `[-0.0472, +0.0760]`
- p-value: `0.6466`
- Status: not met

The fixed LaTeX does not claim that the model beats, outperforms, or is superior to PAM50-ROR.

## Bibliography Checks

`paper/main.bib` includes DOI fields for all 28 references. Each entry also includes a LaTeX `note` with the DOI URL so standard BibTeX output prints a visible DOI link in `paper/main.pdf`. Crossref metadata lookup returned HTTP 200 for all 28 DOI records on 2026-06-02. Key DOI checks:

- Parker 2009 PAM50: `10.1200/JCO.2008.18.1370`
- Perou 2000 molecular portraits: `10.1038/35021093`
- genefu: `10.1093/bioinformatics/btv693`
- TCGA breast cancer: `10.1038/nature11412`
- SCAN-B initiative: `10.1186/s13073-015-0131-9`
- METABRIC Curtis: `10.1038/nature10983`
- METABRIC Pereira: `10.1038/ncomms11479`
- Cox regression: `10.1111/j.2517-6161.1972.tb00899.x`
- Random Survival Forests: `10.1214/08-AOAS169`
- DeepSurv: `10.1186/s12874-018-0482-1`
- Oncotype DX Paik: `10.1056/NEJMoa041588`
- MammaPrint van 't Veer: `10.1038/415530a`
- DerSimonian-Laird: `10.1016/0197-2456(86)90046-2`
- Pencina NRI/IDI: `10.1002/sim.2929`
- Decision curve analysis: `10.1177/0272989X06295361`

## Remaining Caveats

- `pdflatex` itself is not installed, but a local TeX compile was verified with Tectonic 0.16.9 and BibTeX processing.
- The LaTeX source is a reviewer-facing repository artifact. If the journal copy is a Word document, the Word manuscript remains the submitted source of record.
- The 391-variant rescue analysis is described as transparent exploratory post-primary sensitivity analysis, not as a replacement primary endpoint.

## Final Audit Verdict

The specific LaTeX/repo issues from the prior audit are resolved in the active repository. The repo can now honestly be described as containing manuscript LaTeX, a DOI-backed bibliography, compile-safe Table S1 LaTeX, and traceable central results.
