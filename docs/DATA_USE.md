# Data Use, Licensing, and Redistribution

This document states the licensing and redistribution terms for the
BRCA-PathwayML repository and the four source datasets it analyzes. It is
intended for reviewers and downstream users.

## Code license (MIT) — covers CODE ONLY

The repository code is released under the **MIT License** (see `LICENSE`;
Copyright (c) 2026 Suhaan Thayyil and Eshaan Nidee). The MIT license applies to
the **software only** — the Python/R source, analysis scripts, pipeline code,
and documentation. It does **not** grant any rights over the underlying
patient-level datasets, which remain governed by their original providers'
terms (below).

## What is and is not committed in this repository

- **Raw gene-expression matrices are NOT redistributed here.** They are
  obtained from the original providers using the download/parsing scripts and
  the recorded source URLs and SHA256 checksums (see metadata JSON files in
  `data/processed/`).
- What is committed is limited to **de-identified, open-access clinical fields
  and derived/aggregate features** (the seven pathway scores, harmonized
  survival outcomes, and PAM50-ROR comparator scores). These contain **no
  protected health information (PHI)**: no names, dates of birth, full dates of
  service, geographic identifiers, or other HIPAA direct identifiers — only
  cohort-level sample identifiers, age at diagnosis, ordinal stage, receptor
  status, subtype, and survival time/event.
- Users who wish to reproduce the pipeline from raw expression must download the
  source data themselves under each provider's terms.

## Source datasets and their terms

### TCGA-BRCA (training cohort)
- **Access:** Open-access. Clinical and the analyzed expression are open-tier
  TCGA data (controlled-access genotype data is not used).
- **Source:** GDC Data Portal — https://portal.gdc.cancer.gov
  (clinical export also mirrored via UCSC Xena).
- **Terms:** NIH Genomic Data Sharing Policy; open-access tier carries no
  data-use certification requirement. Cite the TCGA Research Network.

### GSE96058 / SCAN-B (external validation)
- **Source:** NCBI GEO accession GSE96058 —
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058
- **Raw file used:** `GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz`
  (download SHA256 `3c717baf7960e1f1477a72744a399a42fabba113c9e7ecb3d55997574d7e9732`,
  per `docs/COHORT_NOTES.md`).
- **Terms:** Use under NCBI GEO terms and the SCAN-B/source release terms.

### METABRIC (external validation; also recurrence/DFS endpoint)
- **Access:** Distributed via cBioPortal under its data-use terms.
- **Source URL used:** `https://datahub.assets.cbioportal.org/brca_metabric.tar.gz`
  (archive SHA256 `6d4683477d6b37a2d7edbedc0df610f67bc456f99e5e1bef6219f37b633a55f7`,
  per `data/processed/03_metabric.metadata.json`). The legacy protocol URL
  `https://cbioportal-datahub.s3.amazonaws.com/brca_metabric.tar.gz` returned
  HTTP 403 on 2026-05-16.
- **cBioPortal study page:**
  https://www.cbioportal.org/study/summary?id=brca_metabric
- **Terms / required citations:** Use under cBioPortal's terms. Cite:
  - Curtis C, et al. *Nature* 2012 (METABRIC discovery cohort).
  - Pereira B, et al. *Nat Commun* 2016 (METABRIC mutational/expression
    landscape).

### GSE20685 (external validation; OS only)
- **Source:** NCBI GEO accession GSE20685 (platform GPL570) —
  https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20685
- **Raw files used:** `GSE20685_series_matrix.txt.gz`
  (SHA256 `e818a4d5834e20bbcf01de515caf856c394abb36a02b9cccc95be05c22d0a279`)
  and GPL570 annotation
  (SHA256 `d7cd44352127b1e34f3a720ebea86093ef255a38f1612a85a2962b71bde8f394`),
  per `data/processed/04_gse20685.metadata.json` and `docs/COHORT_NOTES.md`.
- **Terms:** Use under NCBI GEO terms. Receptor status is not available in the
  parsed series metadata, so this cohort supports overall-survival validation
  only and does not enter the TNBC primary endpoint.

## De-identification statement

The clinical CSVs and parquet files committed to this repository are
**de-identified** and contain **no PHI**. The committed patient-level fields are
limited to a cohort-scoped sample identifier, age at diagnosis (in years),
ordinal/source stage, receptor status, intrinsic subtype, and survival
time/event. No HIPAA direct identifiers are present. Raw expression matrices and
any controlled-access data are not redistributed; they must be obtained from the
original providers under the terms above.

## MSigDB gene sets

Pathway definitions use MSigDB v2024.1.Hs (Hallmark, Reactome, KEGG MEDICUS, PID,
GO BP). Provenance is recorded in
`data/processed/msigdb_v2024_1_Hs.metadata.json`. MSigDB is subject to the Broad
Institute's MSigDB license terms.
