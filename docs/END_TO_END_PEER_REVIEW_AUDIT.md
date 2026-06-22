# End-to-End Peer Review Audit

Generated: 2026-06-02

Scope: GitHub-facing repository, reviewer-facing manuscript source, compiled PDF,
Word export, core result tables, figures, citation metadata, and reproducibility
checks.

## Verdict

The active repository now contains reviewer-facing manuscript artifacts and the
central submitted values are traceable to repository result files. The audit
found and fixed three issues:

1. Standard BibTeX output was not printing DOI fields in the compiled PDF.
2. The Word export required a reproducible title-block cleanup after Pandoc.
3. Stale surrogate-PAM50 notes remained in `src/baselines.py` and
   `results/reports/clinical_baselines_qc.md`.

After remediation, all checks below passed.

## Manuscript Artifacts

Active files:

- `paper/main.tex`
- `paper/main.bib`
- `paper/main.pdf`
- `paper/main.docx`
- `scripts/build_manuscript_docx.py`

Build commands verified:

```bash
cd paper && tectonic --keep-logs main.tex
cd .. && python3 scripts/build_manuscript_docx.py
soffice --headless --convert-to pdf --outdir build/audit/docx_render paper/main.docx
pdftotext paper/main.pdf build/audit/pdf_text/main_pdf.txt
```

Artifact checks:

- `paper/main.pdf`: 29 pages, compiled successfully.
- `paper/main.docx`: rendered successfully through LibreOffice.
- Rendered DOCX PDF: 21 pages.
- DOCX title page visually spot-checked after scripted rebuild.
- Table-heavy DOCX page visually spot-checked.
- Reference DOCX page visually spot-checked.

LaTeX log checks:

- Undefined references: 0.
- Undefined citations: 0.
- LaTeX errors: 0.
- Emergency stop / fatal error: 0.
- Remaining warnings are layout warnings from long table/pathway labels:
  4 overfull hboxes and 26 underfull hboxes.

## Citation Audit

Static citation checks:

- BibTeX entries: 28.
- Inline cite commands: 27.
- Missing cite keys: 0.
- Unused cite keys: 0.
- DOI fields in `paper/main.bib`: 28.
- LaTeX-printable DOI URL notes in `paper/main.bib`: 28.
- Compiled PDF contains visible DOI URLs.

External DOI metadata check:

- Crossref lookup returned HTTP 200 for all 28 DOI records.
- Some publisher redirect targets returned HTTP 403 to scripted requests, but
  Crossref confirmed the DOI records.

## Figure and Table Audit

All figure paths referenced by `paper/main.tex` exist:

- `figures/fig_cohort_flow.pdf`
- `figures/fig_main_forest.pdf`
- `figures/fig_within_subtype_forest.pdf`
- `figures/fig_calibration_plots.pdf`
- `figures/fig_dca_per_cohort.pdf`
- `figures/fig_stability_heatmap.pdf`

Table S1 checks:

- `results/Table_S1_clinical_baselines.tex` uses escaped LaTeX identifiers.
- Stale note saying official genefu should replace surrogate PAM50 was removed.
- Active note now states official genefu PAM50-ROR is primary and the surrogate
  is retained only for historical comparison.

## Number Traceability Audit

Core consistency script:

```bash
python3 scripts/check_submission_consistency.py
```

Status: passed.

Expanded manuscript-number audit:

- 73 table/prose numeric checks against repository CSV/JSON sources.
- Missing values in `paper/main.tex`: 0.

Covered values include:

- Harmonized `N = 4,532`.
- Patient-characteristics evaluable `N = 4,003`.
- Official PAM50-ROR metric-evaluable `N = 4,002`.
- Patient-characteristics events `n = 1,681`.
- Metric-evaluable events `n = 1,680`.
- METABRIC `n = 1,980` patient-characteristics evaluable and `n = 1,979`
  metric-evaluable.
- TNBC primary endpoint `n = 408`, events `n = 214`.
- TNBC meta delta C-index `+0.0144`, 95% CI `-0.0472` to `+0.0760`,
  `p = 0.6466`.
- External per-cohort C-index values for Gradient-Boosted Survival and
  official PAM50-ROR.
- Five-year ICI values and mean ICI values.
- Stability analysis: 100 subsamples, 4,950 pairwise comparisons, mean
  Spearman rho `0.395`.

## Claim-Language Audit

Searches were run for stale or risky strings across active manuscript and
repository-facing documentation.

Resolved:

- No `TODO` or `FIXME` placeholders in `paper/main.tex`.
- No accidental AI attribution strings in active manuscript artifacts.
- No stale "official genefu reproduction should replace this" text remains in
  active source, Table S1, or clinical baseline QC report.

Allowed context:

- Terms such as "superiority", "exceed", and "outperform" remain where they
  describe the pre-registered hypothesis, prior literature, Cox-baseline
  secondary analysis, or the explicitly negative result. The active conclusion
  remains competitive-not-superior.

## Software Checks

Commands:

```bash
python3 -m ruff check src tests scripts
python3 -m pytest tests/
```

Status:

- Ruff: passed.
- Pytest: 3 passed.

Pytest warning:

- `pytest_asyncio` reports a deprecation warning about unset
  `asyncio_default_fixture_loop_scope`. This does not affect the current tests,
  but can be cleaned in future test configuration.

## Remaining Caveats

- The PDF has minor overfull/underfull hbox warnings from long table and
  pathway labels. These are not compile blockers and do not affect citation or
  numeric traceability.
- The DOCX is a generated Word export for reviewer convenience. The PDF/LaTeX
  source remains the cleaner repository manuscript artifact.
- The locked primary endpoint remains not met. The repository should continue
  to avoid claims that pathway ML is superior to PAM50-ROR.

## Final Audit Status

Passed after remediation. The repository is now substantially stronger for peer
review traceability: manuscript source exists, PDF and DOCX build from source,
citations resolve, DOI metadata is visible, central numbers match repository
outputs, and stale PAM50-surrogate wording has been removed.
