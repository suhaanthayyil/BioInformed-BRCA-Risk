# Paper Agent Results Package

Built: 2026-05-17T09:05:42-04:00
Commit SHA: 13bd5f9689efe1a97748edc604c84ed0d477657d

This package is for a paper-writing agent. It contains the numbers, tables, figures, summaries, and manuscript scaffolds needed to write the BMC Cancer version of the paper. It intentionally does not include raw expression data or model binaries.

## Locked Scientific Result

The pre-registered TNBC primary endpoint was NOT MET.

- Headline model: Gradient Boosted Survival
- Comparator: official genefu PAM50-ROR
- TNBC delta C-index: +0.0144
- 95% CI: [-0.0472, +0.0760]
- p-value: 0.6466
- Required threshold: delta >= +0.0300 and p < 0.05

Do not write that the model is superior to PAM50-ROR. The correct framing is comparable or competitive performance with a transparently negative primary endpoint.

## Most Important Files

- docs/results_bullets.md: exact paper-ready result bullets
- docs/PHASE_THREE_SUMMARY.md: stage-by-stage summary and limitations
- docs/STORY.md: locked BMC Cancer framing
- paper/main.tex: BMC manuscript skeleton with Methods and declarations
- paper/main.bib: bibliography
- manuscript/figures/: main figures as PDF and TIFF
- manuscript/tables/: main and supplementary tables as CSV plus Word table bundle
- manuscript/additional_files/: BMC-style supplementary files
- results/: source CSVs and QC reports
- deliverables/SUBMISSION_README.md: remaining author TODOs

## Human TODOs Before Submission

- Write Background, Results prose, Discussion, and Conclusions.
- Fill Suhaan and Eshaan affiliations, emails, ORCID, and contribution details.
- Fill suggested reviewers.
- Install LaTeX and compile paper/main.tex.
- Keep the failed primary endpoint visible and honest.
