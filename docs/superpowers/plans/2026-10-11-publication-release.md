# Publication Release Implementation Plan

> Packaging update (11 October 2026): the manuscript, Supplementary Information,
> and journal correspondence were removed from the public GitHub Release to avoid
> version-of-record and copyright ambiguity. The release now distributes tagged
> source code and machine-readable results only. Asset-publication steps below
> are retained as the historical execution plan.

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Publish a verifiable `v1.2.1` repository release whose default branch, code, results, manuscript, and Supplementary Information all represent the revised 6,644-patient analysis.

**Architecture:** A repository audit script will turn the editor-facing claims into executable checks over committed results, clinical-stage data, version metadata, and the local supplement. Documentation will map each editorial concern to evidence. Release assets will be built from the verified DOCX files, checksummed, attached to a new immutable GitHub Release, and the revised branch will fast-forward the stale default branch.

**Tech Stack:** Python 3.12, pandas, DuckDB, python-docx, pytest, Git, GitHub CLI, LibreOffice-based DOCX renderer.

**Spec:** `docs/PUBLICATION_RELEASE_SPEC.md`

## Global Constraints

- Release version: `1.2.1`; tag: `v1.2.1`; release date: `2026-10-11`.
- Preserve the existing `v1.2.0` tag and history.
- Publish the revised manuscript and Supplementary Information as DOCX and PDF release assets.
- Do not publish the response letter.
- Preserve the current scientific estimates; this task verifies and packages them but does not rerun model selection.
- The default branch must show the revised analysis, not the superseded 4,532-patient analysis.

## Review Focus

- A reader opening the default GitHub page must see 6,644 patients, DeepSurv, and the +0.020 TNBC result.
- The supplement must contain Tables S18-S21 and both Oncotype DX surrogate dispositions.
- Stage must be recorded in TCGA-BRCA, METABRIC, and GSE20685, and absent in GSE96058.
- The release version must agree across `VERSION`, `CITATION.cff`, README, tag, and release title.
- Every public release asset must have a matching SHA256 entry.

---

### Task 1: Executable Publication Audit

**Files:**
- Create: `scripts/check_publication_readiness.py`
- Create: `tests/test_publication_readiness.py`

**Interfaces:**
- Consumes: committed cohort flow JSON, denominator CSV, Tables S18-S21, the unified cohort database, `README.md`, `VERSION`, and `CITATION.cff`.
- Produces: `audit_repository(repo_root: Path, supplement: Path | None = None) -> list[str]` and a CLI that exits nonzero when any publication claim is unsupported.

- [ ] **Step 1: Write failing tests**

Create tests that require version `1.2.1`, the revised headline values, complete registered-endpoint/deviation rows, correct stage availability by cohort, and Supplementary Tables S18-S21 including the Oncotype DX rows.

- [ ] **Step 2: Run the tests and confirm RED**

Run: `.venv/bin/pytest -q tests/test_publication_readiness.py`

Expected: failure because `VERSION` and `check_publication_readiness.py` do not yet exist.

- [ ] **Step 3: Implement the audit**

Read the committed structured outputs directly. Report each failed invariant as a separate line; do not infer readiness from prose alone. Make supplement validation optional so CI can audit the public repository while the release build supplies the local DOCX path.

- [ ] **Step 4: Run the focused tests and confirm GREEN**

Run: `.venv/bin/pytest -q tests/test_publication_readiness.py`

Expected: all focused tests pass after Task 2 supplies the version and documentation metadata.

### Task 2: Version, Chronology, and Evidence Documentation

**Files:**
- Create: `VERSION`
- Create: `docs/PUBLICATION_PACKAGE.md`
- Modify: `README.md`
- Modify: `CITATION.cff`
- Modify: `docs/REGISTRATION_DEVIATIONS.md`
- Modify: `results/Table_S20_registration_deviations.csv`
- Modify: `_submission/Manuscript_prepublication_revision_v1.2.docx`
- Modify: `_submission/BMC_Supplement_prepublication_revision_v1.2.docx`

**Interfaces:**
- Consumes: the release specification and existing Table S18-S21 sources.
- Produces: version `1.2.1`, explicit post-registration chronology, a reviewer evidence map, and synchronized manuscript/supplement wording.

- [ ] **Step 1: Add version and public package metadata**

Set `VERSION` and `CITATION.cff` to `1.2.1` dated 2026-10-11. Add README links to `v1.2.1` and list exactly which cohorts and analyses were new after registration.

- [ ] **Step 2: Add the evidence map**

Document the source file, verification command, and conclusion for registration chronology, S18/S20 completeness, cohort counts, and clinical-stage coding.

- [ ] **Step 3: Synchronize the manuscript and supplement**

Add one explicit chronology sentence to the manuscript Registration paragraph and Supplementary Section S3. Preserve the existing formatting and all scientific values. Fix the duplicated period in Table S20's cohort reason in both structured sources.

- [ ] **Step 4: Render and inspect both DOCX files**

Run the packaged `render_docx.py` for the manuscript and supplement. Inspect every rendered page and iterate on any clipping, overlap, or pagination defect.

- [ ] **Step 5: Run publication and manuscript checks**

Run:

```bash
.venv/bin/pytest -q tests/test_publication_readiness.py
.venv/bin/python scripts/check_publication_readiness.py --supplement _submission/BMC_Supplement_prepublication_revision_v1.2.docx
/Users/suhaan/.cache/codex-runtimes/codex-primary-runtime/dependencies/python/bin/python3 _submission/check_manuscript.py _submission/Manuscript_prepublication_revision_v1.2.docx _submission/BMC_Supplement_prepublication_revision_v1.2.docx
```

Expected: all tests pass, the audit prints no issues, and the manuscript checker reports `0 problem(s)`.

### Task 3: Release Assets and Repository Commit

**Files:**
- Modify: `scripts/calibration_analysis.py`
- Modify: `scripts/dca_analysis.py`
- Modify: `scripts/within_subtype_analysis.py`
- Create: `tests/test_figure_presentation.py`
- Create in a temporary release directory: manuscript DOCX/PDF, supplement DOCX/PDF, `SHA256SUMS.txt`, and `RELEASE_NOTES.md`.

**Interfaces:**
- Consumes: verified DOCX deliverables and the repository at version `1.2.1`.
- Produces: a commit containing version/docs/audit/figure fixes and a checksummed six-file release bundle.

- [ ] **Step 1: Run the full repository gates**

Run: `.venv/bin/pytest -q && .venv/bin/ruff check . && .venv/bin/python scripts/check_submission_consistency.py`

Expected: all tests, lint, and all scientific consistency checks pass.

- [ ] **Step 2: Build release assets**

Render the manuscript and supplement to PDF, copy the source DOCX files under `BRCA-PathwayML_*_v1.2.1` names, write release notes that enumerate the four editorial concerns, and generate SHA256 checksums for every public asset.

- [ ] **Step 3: Commit only intentional repository changes**

Stage the version, README/docs, audit/test, Table S20 typo fix, figure source changes, the three corrected primary figure outputs, and the figure regression test. Exclude the local submission directory, conversation log, lockfile, and incidental regenerated Kaplan-Meier PDFs.

- [ ] **Step 4: Verify the staged commit**

Run the full gates again and inspect `git diff --cached --check` plus the complete staged file list before committing.

### Task 4: Publish and Verify GitHub v1.2.1

**Files:**
- External state: `origin/revision-1-reproducibility`, `origin/main`, tag `v1.2.1`, and GitHub Release `v1.2.1`.

**Interfaces:**
- Consumes: the verified commit and checksummed release bundle.
- Produces: a default branch and immutable GitHub Release that expose the same revised analysis and assets.

- [ ] **Step 1: Push the revision branch and fast-forward main**

Push the current branch, then push it to `main`. Verify the remote commit IDs match and that the update is fast-forward only.

- [ ] **Step 2: Tag and publish the release**

Create annotated tag `v1.2.1`, push it, and use `gh release create` with all manuscript, supplement, checksum, and release-note assets.

- [ ] **Step 3: Wait for GitHub Actions**

Watch the workflow run for the published commit until lint, tests, and consistency complete successfully.

- [ ] **Step 4: Verify live GitHub state**

Use the GitHub API and public pages to confirm: default branch `main` points to the release commit, version `1.2.1` is visible, the stale 4,532-patient summary is gone, every release asset is downloadable, and asset SHA256 values match locally.
