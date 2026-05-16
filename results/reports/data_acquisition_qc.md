# Data Acquisition QC

Timestamp: 2026-05-16T11:16:00-04:00

## Completed Before STOP GATE 1

- GSE96058 FTP listing verified.
- GSE96058 expression gzip downloaded, `gzip -t` passed, and SHA256 recorded.
- GSE96058 expression parsed to parquet with 30,865 gene rows and 3,409 sample columns.
- METABRIC legacy S3 URL returned HTTP 403; current official cBioPortal asset URL was verified with HTTP 200 and downloaded.
- METABRIC archive verified as gzip/tar, tar listing passed, and expected current DataHub files were extracted.
- METABRIC clinical and expression data parsed to parquet.

## Files Generated

- `data/processed/02_gse96058_expression.metadata.json`
- `data/processed/02_gse96058_expression.parquet`
- `data/processed/03_metabric.metadata.json`
- `data/processed/03_metabric_clinical.parquet`
- `data/processed/03_metabric_expression.parquet`

## QC Notes

- Generated parquet files are intentionally ignored by git because they are large derived data.
- Metadata JSON files are retained for provenance and checksums.
- Disk free space after these downloads was approximately 54 GiB, below the original 60 GiB preflight threshold. The threshold was satisfied at Phase 0 start, but future phases may need cleanup or additional disk space.

