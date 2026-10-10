# Raw Data

Raw data files are not stored in this repository due to their size. Use the instructions below to obtain them.

## GSE96058 / SCAN-B Gene Expression

**Source:** [GEO Accession GSE96058](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058)

Download automatically:
```bash
bash download_data.sh
```

Or manually download from:
```
https://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96058/suppl/GSE96058_gene_expression_3273_samples_and_136_replicates_transformed.csv.gz
```

- Compressed: ~565 MB
- Decompressed: ~1.8 GB
- Format: Genes as rows, samples (F1, F2, ...) as columns, log2-transformed expression values

## TCGA-BRCA

**Source:** UCSC Xena TCGA hub (`TCGA.BRCA.sampleMap/HiSeqV2`, `BRCA_clinicalMatrix`) and TCGA-CDR survival (`survival/BRCA_survival.txt`), downloaded by `download_data.sh` into `data/raw/tcga/`.

- Expression: log2(RSEM normalized count + 1), 20,530 genes
- Parsed by `scripts/ingest_raw.py` into `data/processed/01_tcga_expression.parquet` and `data/clinical/tcga_brca_clinical.csv`

## GSE96058 clinical

The two GEO series matrices (`GSE96058-GPL11154`, `GSE96058-GPL18573`) are downloaded into `data/raw/gse96058/` and parsed by `scripts/ingest_raw.py` into `data/clinical/gse96058_clinical.csv` (3,273 patients, replicates removed).

## Clinical Data

The legacy v1 (preprint) cohort files are in `data/legacy_v1_preprint/` (see its README).
